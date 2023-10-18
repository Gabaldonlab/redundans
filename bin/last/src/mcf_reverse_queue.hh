// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef MCF_REVERSE_QUEUE
#define MCF_REVERSE_QUEUE

#include <new>
#include <stdlib.h>
#include <string.h>

namespace mcf {

template <typename T> class ReverseQueue {
public:
  ReverseQueue() : buf(0), beg(0), end(0) {}

  void clear() {
    beg = end;
  }

  void push(const T &item, unsigned numOfOldItemsToKeep) {
    if (beg == buf) {
      size_t len = end - beg;
      if (numOfOldItemsToKeep < len / 2) {
	beg = end - numOfOldItemsToKeep;
	memcpy(beg, buf, numOfOldItemsToKeep * sizeof(T));
      } else {
	size_t newLen = len * 2 + 128;
	size_t newBytes = newLen * sizeof(T);
	if (newBytes <= len * sizeof(T)) throw std::bad_alloc();
	void *p = malloc(newBytes);
	if (!p) throw std::bad_alloc();
	end = static_cast<T *>(p) + newLen;
	beg = end - numOfOldItemsToKeep;
	if (buf) memcpy(beg, buf, numOfOldItemsToKeep * sizeof(T));
	free(buf);
	buf = static_cast<T *>(p);
      }
    }
    --beg;
    *beg = item;
  }

  const T *begin() const {
    return beg;
  }

  ~ReverseQueue() {
    free(buf);
  }

private:
  T *buf;
  T *beg;
  T *end;
};

}

#endif
