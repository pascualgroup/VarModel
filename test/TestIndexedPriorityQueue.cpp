//
//  TestIndexedPriorityQueue.cpp
//  malariamodel
//
//  Created by Ed Baskerville on 2/21/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#include "catch.hpp"
#include "types.h"
#include "IndexedPriorityQueue.h"

using namespace std;

static int compare(size_t i1, size_t i2)
{
	return i1 < i2 ? -1 : (i1 > i2 ? 1 : 0);
}

TEST_CASE("IndexedPriorityQueue behaves properly", "[IndexedPriorityQueue]")
{
	rng_t rng(100);
	
	SECTION("Simple initial heapify.") {
		size_t preHeap[] = {6, 0, 1, 2, 3, 4, 5};
		IndexedPriorityQueueProbe probe(vector<size_t>(begin(preHeap), end(preHeap)), compare);
		probe.heapifyDown(0);
		
		size_t postHeap[] = {0, 2, 1, 6, 3, 4, 5};
		REQUIRE(probe.getHeap() == vector<size_t>(begin(postHeap), end(postHeap)));
	}
	
	SECTION("Incomplete heapify.") {
		size_t preHeap[] = {3, 0, 1, 2};
		IndexedPriorityQueueProbe probe(vector<size_t>(begin(preHeap), end(preHeap)), compare);
		probe.heapifyDown(0);
		
		size_t postHeap[] = {0, 2, 1, 3};
		REQUIRE(probe.getHeap() == vector<size_t>(begin(postHeap), end(postHeap)));
	}
	
	SECTION("Nonzero index with nothing to do.") {
		size_t preHeap[] = {5, 0, 1, 2, 3, 4};
		IndexedPriorityQueueProbe probe(vector<size_t>(begin(preHeap), end(preHeap)), compare);
		probe.heapifyDown(3);
		
		size_t postHeap[] = {5, 0, 1, 2, 3, 4};
		REQUIRE(probe.getHeap() == vector<size_t>(begin(postHeap), end(postHeap)));
	}
	
	SECTION("Nonzero index with something to do.") {
		size_t preHeap[] = {0, 6, 1, 2, 3, 4, 5};
		IndexedPriorityQueueProbe probe(vector<size_t>(begin(preHeap), end(preHeap)), compare);
		probe.heapifyDown(1);
		
		size_t postHeap[] = {0, 2, 1, 6, 3, 4, 5};
		REQUIRE(probe.getHeap() == vector<size_t>(begin(postHeap), end(postHeap)));
	}
	
	SECTION("Borderline index heapify down.") {
		size_t preHeap[] = {0, 1, 6, 2, 3, 4, 5};
		IndexedPriorityQueueProbe probe(vector<size_t>(begin(preHeap), end(preHeap)), compare);
		probe.heapifyDown(2);
		
		size_t postHeap[] = {0, 1, 4, 2, 3, 6, 5};
		REQUIRE(probe.getHeap() == vector<size_t>(begin(postHeap), end(postHeap)));
	}
	
	SECTION("Heapify up.") {
		size_t preHeap[] = {1, 2, 3, 4, 0, 5, 6};
		IndexedPriorityQueueProbe probe(vector<size_t>(begin(preHeap), end(preHeap)), compare);
		probe.heapifyUp(4);
		
		size_t postHeap[] = {0, 1, 3, 4, 2, 5, 6};
		REQUIRE(probe.getHeap() == vector<size_t>(begin(postHeap), end(postHeap)));
	}
	
	SECTION("Reverse heap build.") {
		size_t preHeap[] = {7, 6, 5, 4, 3, 2, 1, 0};
		IndexedPriorityQueueProbe probe(vector<size_t>(begin(preHeap), end(preHeap)), compare);
		probe.buildHeap();
		
		size_t postHeap[] = {0, 3, 1, 4, 7, 2, 5, 6};
		REQUIRE(probe.getHeap() == vector<size_t>(begin(postHeap), end(postHeap)));
	}
	
	SECTION("Update up.") {
		int values[] = {1, 2, 3, 4, 5, 6, 7};
		
		size_t preHeap[] = {0, 1, 2, 3, 4, 5, 6};
		
		IndexedPriorityQueueProbe probe(vector<size_t>(begin(preHeap), end(preHeap)),
			[&values](size_t i1, size_t i2) {
				return values[i1] < values[i2] ? -1 : (values[i1] > values[i2] ? 1 : 0);
			}
		);
		values[5] = 0;
		probe.update(5);
		
		size_t postHeap[] = {5, 1, 0, 3, 4, 2, 6};
		REQUIRE(probe.getHeap() == vector<size_t>(begin(postHeap), end(postHeap)));
	}
	
	SECTION("Update down.") {
		int values[] = {1, 2, 3, 4, 6, 7, 8};
		
		size_t preHeap[] = {0, 1, 2, 3, 4, 5, 6};
		
		IndexedPriorityQueueProbe probe(vector<size_t>(begin(preHeap), end(preHeap)),
			[&values](size_t i1, size_t i2) {
				return values[i1] < values[i2] ? -1 : (values[i1] > values[i2] ? 1 : 0);
			}
		);
		values[0] = 5;
		probe.update(0);
		
		size_t postHeap[] = {1, 3, 2, 0, 4, 5, 6};
		REQUIRE(probe.getHeap() == vector<size_t>(begin(postHeap), end(postHeap)));
	}
}
