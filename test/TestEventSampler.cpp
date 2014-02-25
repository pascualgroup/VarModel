//
//  TestEventSampler.cpp
//  malariamodel
//
//  Created by Ed Baskerville on 2/21/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#include "catch.hpp"
#include "random.h"
#include "EventSampler.h"

using namespace std;

TEST_CASE("EventSampler behaves properly", "[EventSampler]")
{
	rng_t rng(100);
	
	SECTION("Approximately the right number of events for rate = 1.") {
		unique_ptr<Event> event(new SimpleEvent(1.0));
		vector<Event *> events;
		events.push_back(event.get());
		
		EventSampler sampler(events, rng);
		while(sampler.getTime() < 1000.0)
		{
			Event * event;
			size_t eventId;
			double dt;
			sampler.performNextEvent(event, eventId, dt);
		}
		REQUIRE(sampler.getEventCount() > 900);
		REQUIRE(sampler.getEventCount() < 1100);
		REQUIRE(sampler.getEventCount() == event->getCount());
	}
	
	SECTION("Approximately the right number of events for rate = 2.") {
		unique_ptr<Event> event(new SimpleEvent(2.0));
		vector<Event *> events;
		events.push_back(event.get());
		
		EventSampler sampler(events, rng);
		while(sampler.getTime() < 500.0)
		{
			Event * event;
			size_t eventId;
			double dt;
			sampler.performNextEvent(event, eventId, dt);
		}
		REQUIRE(sampler.getEventCount() > 900);
		REQUIRE(sampler.getEventCount() < 1100);
		REQUIRE(sampler.getEventCount() == event->getCount());
	}
	
	SECTION("Approximately the right number of events for two equal-rate events") {
		vector<Event *> events;
		unique_ptr<Event> event1(new SimpleEvent(1.0));
		events.push_back(event1.get());
		unique_ptr<Event> event2(new SimpleEvent(1.0));
		events.push_back(event2.get());
		
		EventSampler sampler(events, rng);
		while(sampler.getTime() < 1000.0)
		{
			Event * event;
			size_t eventId;
			double dt;
			sampler.performNextEvent(event, eventId, dt);
		}
		REQUIRE(event1->getCount() > 900);
		REQUIRE(event1->getCount() < 1100);
		REQUIRE(event2->getCount() > 900);
		REQUIRE(event2->getCount() < 1100);
		REQUIRE(sampler.getEventCount() == event1->getCount() + event2->getCount());
	}
	
	SECTION("Approximately the right number of events for two different-rate events") {
		vector<Event *> events;
		unique_ptr<Event> event1(new SimpleEvent(2.0));
		events.push_back(event1.get());
		unique_ptr<Event> event2(new SimpleEvent(1.0));
		events.push_back(event2.get());
		
		EventSampler sampler(events, rng);
		while(sampler.getTime() < 1000.0)
		{
			Event * event;
			size_t eventId;
			double dt;
			sampler.performNextEvent(event, eventId, dt);
		}
		REQUIRE(event1->getCount() > 1800);
		REQUIRE(event1->getCount() < 2200);
		REQUIRE(event2->getCount() > 900);
		REQUIRE(event2->getCount() < 1100);
		REQUIRE(sampler.getEventCount() == event1->getCount() + event2->getCount());
	}
	
	SECTION("Approximately the right number of events before and after rate change") {
		SimpleEvent event(1.0);
		vector<Event *> events(1, &event);
		
		EventSampler sampler(events, rng);
		while(sampler.getTime() < 1000.0)
		{
			Event * event;
			size_t eventId;
			double dt;
			sampler.performNextEvent(event, eventId, dt);
		}
		REQUIRE(event.getCount() > 900);
		REQUIRE(event.getCount() < 1100);
		
		sampler.resetEventCounts();
		
		event.setRate(2.0);
		sampler.updateEvent(&event);
		while(sampler.getTime() < 1500.0)
		{
			Event * event;
			size_t eventId;
			double dt;
			sampler.performNextEvent(event, eventId, dt);
		}
		REQUIRE(event.getCount() > 900);
		REQUIRE(event.getCount() < 1100);
		REQUIRE(sampler.getEventCount() == event.getCount());
	}
	
	SECTION("Approximately the right number of events before and after several rate changes with two events.") {
		SimpleEvent event1(1.0);
		SimpleEvent event2(1.0);
		vector<Event *> events;
		events.push_back(&event1);
		events.push_back(&event2);
		
		EventSampler sampler(events, rng);
		
		Event * event;
		size_t eventId;
		double dt;
		while(sampler.getTime() < 1000.0)
		{
			sampler.performNextEvent(event, eventId, dt);
		}
		REQUIRE(event1.getCount() > 900);
		REQUIRE(event1.getCount() < 1100);
		REQUIRE(event2.getCount() > 900);
		REQUIRE(event2.getCount() < 1100);
		REQUIRE(sampler.getEventCount() == event1.getCount() + event2.getCount());
		
		sampler.resetEventCounts();
		event1.setRate(2.0);
		sampler.updateEvent(&event1);
		
		while(sampler.getTime() < 2000.0)
		{
			sampler.performNextEvent(event, eventId, dt);
		}
		REQUIRE(event1.getCount() > 1800);
		REQUIRE(event1.getCount() < 2200);
		REQUIRE(event2.getCount() > 900);
		REQUIRE(event2.getCount() < 1100);
		
		sampler.resetEventCounts();
		event2.setRate(2.0);
		sampler.updateEvent(&event2);
		
		while(sampler.getTime() < 2500.0)
		{
			sampler.performNextEvent(event, eventId, dt);
		}
		REQUIRE(event1.getCount() > 900);
		REQUIRE(event1.getCount() < 1100);
		REQUIRE(event2.getCount() > 900);
		REQUIRE(event2.getCount() < 1100);
		
		sampler.resetEventCounts();
		event1.setRate(0.5);
		sampler.updateEvent(&event1);
		event2.setRate(0.5);
		sampler.updateEvent(&event2);
		
		while(sampler.getTime() < 4500.0)
		{
			sampler.performNextEvent(event, eventId, dt);
		}
		REQUIRE(event1.getCount() > 900);
		REQUIRE(event1.getCount() < 1100);
		REQUIRE(event2.getCount() > 900);
		REQUIRE(event2.getCount() < 1100);
		
		sampler.resetEventCounts();
		event1.setRate(0.5);
		sampler.updateEvent(&event1);
		event2.setRate(5.0);
		sampler.updateEvent(&event2);
		
		while(sampler.getTime() < 6500.0)
		{
			sampler.performNextEvent(event, eventId, dt);
		}
		REQUIRE(event1.getCount() > 900);
		REQUIRE(event1.getCount() < 1100);
		REQUIRE(event2.getCount() > 9000);
		REQUIRE(event2.getCount() < 11000);
	}
	
	SECTION("Approximately the right number of events before and after many rate changes with many events.") {
		vector<SimpleEvent> events(101, SimpleEvent(1.0));
		vector<Event *> eventPtrs;
		for(int i = 0; i < events.size(); i++) {
			eventPtrs.push_back(&events[i]);
		}
		
		EventSampler sampler(eventPtrs, rng);
		uniform_int_distribution<size_t> ud(0, 100);
		size_t specialId = 0;
		for(int i = 0; i < 100; i++) {
			events[specialId].setRate(1.0);
			sampler.updateEvent(eventPtrs[specialId]);
			specialId = ud(rng);
			events[specialId].setRate(100.0);
			sampler.updateEvent(eventPtrs[specialId]);
			while(sampler.getTime() < 10.0 * (i+1)) {
				Event * event;
				size_t eventId;
				double dt;
				sampler.performNextEvent(event, eventId, dt);
			}
			REQUIRE(sampler.getEventCount() > 1800);
			REQUIRE(sampler.getEventCount() < 2200);
			REQUIRE(events[specialId].getCount() > 900);
			REQUIRE(events[specialId].getCount() < 1100);
			sampler.resetEventCounts();
		}
	}
}
