// license:GPLv3+
//
// A part's values by their lifetime (runs/2026-11-08-state-split/design.md): the part declares them, this shared code alone saves, restores and wipes them.
//   state     a value a later step, cycle or retake can read before it is written again: a plain struct the part derives from, copied whole by the engine's snapshot and the net's retake
//   vectors   a state that is a list sized per island (the cord's touches): its size and items, copied whole
//   queues    events kept until the owner reads them: their counts are saved, a restore cuts each back to its count (the replay makes the later ones again)
//   step work a value written and read inside one Group::Step (all its attempts): wiped to its defaults when the step starts
//   attempt   a value written and read inside one attempt (BeginStep .. Eval .. Accept): wiped to its defaults at the attempt's start
//   use work  the group step's own: everything it writes while it steps one island (or one part alone), wiped at the start of the next use (Group::BeginUse)
// A part's settings and its shape are none of these: they are not saved. Records are kept by their offset from this object, so a copied part keeps its own.

#pragma once

#include <cassert>
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <vector>

// A work record back to its defaults (one rule for every kind of record): a plain struct or value by its default construction; a list emptied with its storage kept;
// a list of lists each emptied. A kind of its own (the group's buffers) adds its own overload beside it
template <class W> void WipeRecord(W &w) { w = W(); }
template <class T> void WipeRecord(std::vector<T> &v) { v.clear(); }
template <class T> void WipeRecord(std::vector<std::vector<T>> &v)
{
   for (std::vector<T> &l : v)
      l.clear();
}

// A list of at most N trivially copyable items held in place: a part's state can hold it (a body's stops)
template <class T, int N> struct FixedList
{
   static_assert(std::is_trivially_copyable<T>::value, "FixedList holds plain items");
   T item[N];
   int n = 0;
   int size() const { return n; }
   bool empty() const { return n == 0; }
   void clear() { n = 0; }
   void push_back(const T &v)
   {
      assert(n < N);
      item[n++] = v;
   }
   T &operator[](const int i) { return item[i]; }
   const T &operator[](const int i) const { return item[i]; }
   T *begin() { return item; }
   T *end() { return item + n; }
   const T *begin() const { return item; }
   const T *end() const { return item + n; }
};

class PartRecords
{
public:
   // ---- what a part declares, once, in its constructor ----
   template <class S> void KeepState(S *s)
   {
      static_assert(std::is_trivially_copyable<S>::value, "a state record is plain bytes: no pointer to its own storage, no vector");
      assert(nState < kMax);
      state[nState++] = { Off(s), sizeof(S), nullptr };
   }
   template <class T> void KeepVector(std::vector<T> *v)
   {
      static_assert(std::is_trivially_copyable<T>::value, "a kept list holds plain items");
      assert(nVec < kMax);
      vec[nVec++] = { Off(v), sizeof(T), [](const void *p) { return static_cast<const std::vector<T> *>(p)->size(); },
         [](void *p, const size_t n) { static_cast<std::vector<T> *>(p)->resize(n); return (void *)static_cast<std::vector<T> *>(p)->data(); },
         [](const void *p) { return (const void *)static_cast<const std::vector<T> *>(p)->data(); } };
   }
   template <class Q> void KeepQueue(std::vector<Q> *q)
   {
      assert(nQueue < kMax);
      queue[nQueue++] = { Off(q), [](const void *p) { return static_cast<const std::vector<Q> *>(p)->size(); },
         [](void *p, const size_t n) { std::vector<Q> &v = *static_cast<std::vector<Q> *>(p); if (v.size() > n) v.resize(n); } };
   }
   template <class W> void WipeEachStep(W *w)
   {
      assert(nStep < kMax);
      step[nStep++] = { Off(w), sizeof(W), [](void *p) { WipeRecord(*static_cast<W *>(p)); } };
   }
   template <class W> void WipeEachAttempt(W *w)
   {
      assert(nAttempt < kMax);
      attempt[nAttempt++] = { Off(w), sizeof(W), [](void *p) { WipeRecord(*static_cast<W *>(p)); } };
   }
   template <class W> void WipeEachUse(W *w) { use.push_back({ Off(w), sizeof(W), [](void *p) { WipeRecord(*static_cast<W *>(p)); } }); }

   // ---- the shared code: the snapshot, the net's retake, the step's start ----
   size_t StateBytes() const
   {
      size_t n = 0;
      for (int i = 0; i < nState; i++)
         n += state[i].size;
      for (int i = 0; i < nVec; i++)
         n += sizeof(size_t) + vec[i].count(At(vec[i].off)) * vec[i].item;
      return n + nQueue * sizeof(size_t);
   }
   // the state records, the kept lists, then the queues' counts, at out (StateBytes of them)
   void SaveRecords(unsigned char *out) const
   {
      for (int i = 0; i < nState; i++)
      {
         std::memcpy(out, At(state[i].off), state[i].size);
         out += state[i].size;
      }
      for (int i = 0; i < nVec; i++)
      {
         const size_t c = vec[i].count(At(vec[i].off)), bytes = c * vec[i].item;
         std::memcpy(out, &c, sizeof c);
         out += sizeof c;
         if (bytes > 0)
            std::memcpy(out, vec[i].data(At(vec[i].off)), bytes);
         out += bytes;
      }
      for (int i = 0; i < nQueue; i++)
      {
         const size_t c = queue[i].count(At(queue[i].off));
         std::memcpy(out, &c, sizeof c);
         out += sizeof c;
      }
   }
   void RestoreRecords(const unsigned char *in)
   {
      for (int i = 0; i < nState; i++)
      {
         std::memcpy(At(state[i].off), in, state[i].size);
         in += state[i].size;
      }
      for (int i = 0; i < nVec; i++)
      {
         size_t c;
         std::memcpy(&c, in, sizeof c);
         in += sizeof c;
         void *const d = vec[i].resize(At(vec[i].off), c);
         if (c > 0)
            std::memcpy(d, in, c * vec[i].item);
         in += c * vec[i].item;
      }
      for (int i = 0; i < nQueue; i++)
      {
         size_t c;
         std::memcpy(&c, in, sizeof c);
         in += sizeof c;
         queue[i].cut(At(queue[i].off), c);
      }
   }
   void WipeStepWork()
   {
      for (int i = 0; i < nStep; i++)
         step[i].wipe(At(step[i].off));
   }
   void WipeAttemptWork()
   {
      for (int i = 0; i < nAttempt; i++)
         attempt[i].wipe(At(attempt[i].off));
   }
   void WipeUseWork()
   {
      for (const Rec &r : use)
         r.wipe(At(r.off));
   }

private:
   static constexpr int kMax = 4;
   struct Rec { ptrdiff_t off; size_t size; void (*wipe)(void *); };
   struct Queue { ptrdiff_t off; size_t (*count)(const void *); void (*cut)(void *, size_t); };
   struct Vec { ptrdiff_t off; size_t item; size_t (*count)(const void *); void *(*resize)(void *, size_t); const void *(*data)(const void *); };
   Rec state[kMax] = {}, step[kMax] = {}, attempt[kMax] = {};
   Queue queue[kMax] = {};
   Vec vec[kMax] = {};
   std::vector<Rec> use;   // (the group's many buffers)
   int nState = 0, nStep = 0, nAttempt = 0, nQueue = 0, nVec = 0;
   ptrdiff_t Off(const void *p) const { return static_cast<const char *>(p) - reinterpret_cast<const char *>(this); }
   void *At(const ptrdiff_t off) { return reinterpret_cast<char *>(this) + off; }
   const void *At(const ptrdiff_t off) const { return reinterpret_cast<const char *>(this) + off; }
};
