#ifndef __malariamodel__IndexedPriorityQueue__
#define __malariamodel__IndexedPriorityQueue__

#include <vector>
#include <functional>

class IndexedPriorityQueueProbe;

class IndexedPriorityQueue
{
friend class IndexedPriorityQueueProbe;
public:
	IndexedPriorityQueue(size_t size, std::function<int(size_t o1, size_t o2)> compare);
	void update(size_t id);
	size_t getHead();
private:
	IndexedPriorityQueue(std::vector<size_t> initialHeap, std::function<int(size_t o1, size_t o2)> compare);
	std::function<int(size_t o1, size_t o2)> compare;
	std::vector<size_t> heap;
	std::vector<size_t> indexes;
	
	std::vector<size_t> makeIndexes();
	
	void buildHeap();
	bool heapifyUp(size_t i);
	bool heapifyDown(size_t i);
	void swap(size_t i1, size_t i2);
};

class IndexedPriorityQueueProbe
{
public:
	IndexedPriorityQueueProbe(size_t size, std::function<int(size_t o1, size_t o2)> compare) :
		queue(size, compare)
	{
	}
	
	IndexedPriorityQueueProbe(std::vector<size_t> initialHeap, std::function<int(size_t o1, size_t o2)> compare) :
		queue(initialHeap, compare)
	{
	}
	
	IndexedPriorityQueue queue;
	
	void update(size_t id) { queue.update(id); }
	size_t getHead() { return queue.getHead(); }
	
	void buildHeap() { queue.buildHeap(); }
	bool heapifyUp(size_t i) { return queue.heapifyUp(i); }
	bool heapifyDown(size_t i) { return queue.heapifyDown(i); }
	void swap(size_t i1, size_t i2) { queue.swap(i1, i2); }
	
	std::vector<size_t> getHeap() { return queue.heap; }
};

#endif
