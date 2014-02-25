#ifndef __malariamodel__EventSampler__
#define __malariamodel__EventSampler__

#include <vector>
#include "random.h"
#include "IndexedPriorityQueue.h"

class Event;

typedef std::function<std::vector<Event *>(double)> PerformEventFunction;
typedef std::function<double()> GetRateFunction;

class Event
{
friend class EventSampler;
public:
	Event() : count(0) {}
	
	virtual double getRate() = 0;
	virtual std::vector<Event *> performEvent(double time) = 0;
	
	uint64_t getCount() { return count; }
	
private:
	size_t id;
	double rate;
	double time;
	uint64_t count;
};

class SimpleEvent : public Event
{
public:
	SimpleEvent(SimpleEvent const & event) : Event(), constantRate(event.constantRate) {}
	
	SimpleEvent(double rate) : Event(), constantRate(rate) {}
	virtual void setRate(double rate) { constantRate = rate; }
	virtual double getRate() { return constantRate; }
	virtual std::vector<Event *> performEvent(double time) { return std::vector<Event *>(0); }
private:
	double constantRate;
};

class FunctionEvent : public Event
{
friend class EventSampler;
public:
	FunctionEvent() : Event(), performEventFunc(NULL), getRateFunc(NULL) {}
	
	FunctionEvent(PerformEventFunction performEventFunc, GetRateFunction getRateFunc)
	: Event(), performEventFunc(performEventFunc), getRateFunc(getRateFunc) {}
	
	virtual double getRate()
	{
		return getRateFunc();
	}
	
	virtual std::vector<Event *> performEvent(double time)
	{
		if(performEventFunc == NULL) {
			return std::vector<Event *>(0);
		}
		return performEventFunc(time);
	}
	
	PerformEventFunction performEventFunc;
	GetRateFunction getRateFunc;
};

class EventSampler
{
public:
	EventSampler(std::vector<Event *> const & events, rng_t & rng);
	void updateEvent(Event * event);
	
	bool performNextEvent(Event * & eventOut, size_t & eventIdOut, double & dtOut);
	double getTime();
	
	uint64_t getEventCount();
	void resetEventCounts();

private:
	rng_t * rngPtr;
	
	double time;
	uint64_t eventCount;
	
	// Events
	std::vector<Event *> events;
	
	// Indexed priority queue of events sorted by times
	IndexedPriorityQueue queue;
	
	std::vector<Event *> initializeEvents(std::vector<Event *> const & events);
	
	// Comparison function used by queue
	int queueCompare(size_t e1, size_t e2);
};

#endif
