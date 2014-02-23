#include "EventSampler.h"

#include <vector>
#include <random>
#include <cassert>

using namespace std;

EventSampler::EventSampler(std::vector<Event *> const & events, rng_t & rng)
:
	rngPtr(&rng),
	time(0.0),
	eventCount(0),
	events(initializeEvents(events)),
	queue(events.size(), std::bind(
		&EventSampler::queueCompare, this, placeholders::_1, placeholders::_2
	))
{
}

std::vector<Event *> EventSampler::initializeEvents(std::vector<Event *> const & events)
{
	for(size_t i = 0; i < events.size(); i++)
	{
		events[i]->id = i;
		events[i]->rate = events[i]->getRate();
		exponential_distribution<> expDist(events[i]->rate);
		events[i]->time = expDist(*rngPtr);
	}
	return events;
}

int EventSampler::queueCompare(size_t e1, size_t e2)
{
	return events[e1]->time < events[e2]->time ? -1 : (events[e1]->time > events[e2]->time ? 1 : 0);
}

void EventSampler::updateEvent(Event * event)
{
	double rate = event->getRate();
	assert(rate >= 0);
	double oldRate = event->rate;
	
	if(oldRate != rate) {
		double oldTime = event->time;
		event->rate = rate;
		
		// If the new rate is 0, then the event will never happen.
		if(rate == 0) {
			event->time = std::numeric_limits<double>::infinity();
		}
		// If the old rate was 0, then we need to sample a new time from scratch.
		else if(oldRate == 0) {
			exponential_distribution<> timeDist(rate);
			event->time = time + timeDist(*rngPtr);
		}
		// If the old rate was nonzero, we can rescale the previously generated time.
		else {
			event->time = time + (oldTime - time) * oldRate / rate;
		}
		queue.update(event->id);
	}
}

bool EventSampler::performNextEvent(Event * & eventOut, size_t & eventIdOut, double & dtOut)
{
	size_t eventId = queue.getHead();
	Event * event = events[eventId];
	dtOut = event->time - time;
	time = event->time;
	
	// If all event times are infinite, the simulation is done: return false.
	if(isinf(time)) {
		return false;
	}
	// Otherwise get the next event and resample its time.
	else {
		eventIdOut = eventId;
		eventOut = event;
		exponential_distribution<> timeDist(event->rate);
		double dtFuture = timeDist(*rngPtr);
		event->time = time + dtFuture;
		queue.update(eventId);
		
		vector<Event *> changedEvents = event->performEvent(time);
		for(Event * changedEvent : changedEvents) {
			updateEvent(changedEvent);
		}
		
		eventCount++;
		event->count++;
		
		return true;
	}
}

double EventSampler::getTime()
{
	return time;
}

uint64_t EventSampler::getEventCount()
{
	return eventCount;
}

void EventSampler::resetEventCounts()
{
	eventCount = 0;
	for(Event * event : events) {
		event->count = 0;
	}
}
