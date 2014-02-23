#include "IndexedPriorityQueue.h"

#include <cassert>
#include "vecutils.h"

using namespace std;

static size_t getParent(size_t i)
{
	assert(i > 0);
	return (i-1)/2;
}

static size_t getLeft(size_t i)
{
	return 2*i + 1;
}

IndexedPriorityQueue::IndexedPriorityQueue(
	size_t size,
	std::function<int(size_t o1, size_t o2)> compare
):
	compare(compare),
	heap(makeRange(size)),
	indexes(makeIndexes())
{
	assert(size > 0);
	
	buildHeap();
}

IndexedPriorityQueue::IndexedPriorityQueue(
	std::vector<size_t> initialHeap,
	std::function<int(size_t o1, size_t o2)> compare
):
	compare(compare),
	heap(initialHeap),
	indexes(makeIndexes())
{
	assert(heap.size() > 0);
}

std::vector<size_t> IndexedPriorityQueue::makeIndexes()
{
	vector<size_t> tmpIndexes(heap.size(), numeric_limits<size_t>::max());
	for(size_t i = 0; i < heap.size(); i++)
	{
		assert(tmpIndexes[heap[i]] == numeric_limits<size_t>::max());
		indexes[heap[i]] = i;
	}
	
	return tmpIndexes;
}

void IndexedPriorityQueue::update(size_t id)
{
	size_t index = indexes[id];
	if(!heapifyUp(index))
	{
		heapifyDown(index);
	}
}

size_t IndexedPriorityQueue::getHead()
{
	return heap[0];
}

void IndexedPriorityQueue::buildHeap()
{
	if(heap.size() < 2) {
		return;
	}
	
	// Build heap
	size_t maxParentIndex = (heap.size() - 2)/2;
	for(size_t i = maxParentIndex + 1; i > 0; i--) {
		heapifyDown(i - 1);
	}
}

bool IndexedPriorityQueue::heapifyUp(size_t i)
{
	assert(i  < heap.size());
	
	bool moved = false;
	while(i > 0) {
		size_t parent = getParent(i);
		if(compare(heap[i], heap[parent]) < 0) {
			swap(i, parent);
			i = parent;
			moved = true;
		}
		else break;
	}
	return moved;
}

bool IndexedPriorityQueue::heapifyDown(size_t i)
{
	assert(i < heap.size());
	if(heap.size() < 2) {
		return false;
	}
	
	size_t nullIndex = std::numeric_limits<size_t>::max();
	
	bool moved = false;
	
	size_t size = heap.size();
	size_t maxParentIndex = (size - 2) / 2;
	while(i <= maxParentIndex)
	{
		size_t nVal = heap[i];
		
		size_t left = getLeft(i);
		if(left >= size) break;
		size_t leftVal = heap[left];
		
		size_t right = left + 1;
		size_t rightVal = right < size ? heap[right] : nullIndex;
		
		size_t min = i;
		
		if(compare(leftVal, nVal) < 0)
		{
			if(rightVal != nullIndex && compare(rightVal, leftVal) < 0) {
				min = right;
			}
			else {
				min = left;
			}
		}
		else if(rightVal != nullIndex && compare(rightVal, nVal) < 0) {
			min = right;
		}
		else {
			break;
		}
		
		moved = true;
		swap(i, min);
		i = min;
	}
	return moved;
}

void IndexedPriorityQueue::swap(size_t i1, size_t i2)
{
	assert(i1 < heap.size());
	assert(i2 < heap.size());
	assert(i1 != i2);
	
	size_t e1 = heap[i2];
	size_t e2 = heap[i1];
	heap[i1] = e1;
	heap[i2] = e2;
	
	indexes[e1] = i1;
	indexes[e2] = i2;
}
