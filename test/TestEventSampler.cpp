//
//  TestEventSampler.cpp
//  multistrain_abm
//
//  Created by Ed Baskerville on 2/21/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#include "catch.hpp"
#include "types.h"
#include "EventSampler.h"

using namespace std;

TEST_CASE("EventSampler behaves properly", "[EventSampler]")
{
	rng_t rng(100);
	
	SECTION("Approximately the right number of events for rate = 1.") {
		EventSampler sampler(vector<double>(1, 1.0), rng);
		while(sampler.getTime() < 1000.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
		}
		REQUIRE(sampler.getEventCount() > 900);
		REQUIRE(sampler.getEventCount() < 1100);
		REQUIRE(sampler.getEventCount() == sampler.getEventCount(0));
	}
	
	SECTION("Approximately the right number of events for rate = 2.") {
		EventSampler sampler(vector<double>(1, 2.0), rng);
		while(sampler.getTime() < 500.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
		}
		REQUIRE(sampler.getEventCount() > 900);
		REQUIRE(sampler.getEventCount() < 1100);
		REQUIRE(sampler.getEventCount() == sampler.getEventCount(0));
	}
	
	SECTION("Approximately the right number of events for two equal-rate events") {
		EventSampler sampler(vector<double>(2, 1.0), rng);
		while(sampler.getTime() < 1000.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
		}
		REQUIRE(sampler.getEventCount(0) > 900);
		REQUIRE(sampler.getEventCount(0) < 1100);
		REQUIRE(sampler.getEventCount(1) > 900);
		REQUIRE(sampler.getEventCount(1) < 1100);
		REQUIRE(sampler.getEventCount() == sampler.getEventCount(0) + sampler.getEventCount(1));
	}
	
	SECTION("Approximately the right number of events for two different-rate events") {
		double initialRates[] = {2.0, 1.0};
		EventSampler sampler(vector<double>(begin(initialRates), end(initialRates)), rng);
		while(sampler.getTime() < 1000.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
		}
		REQUIRE(sampler.getEventCount(0) > 1800);
		REQUIRE(sampler.getEventCount(0) < 2200);
		REQUIRE(sampler.getEventCount(1) > 900);
		REQUIRE(sampler.getEventCount(1) < 1100);
		REQUIRE(sampler.getEventCount() == sampler.getEventCount(0) + sampler.getEventCount(1));
	}
	
	SECTION("Approximately the right number of events before and after rate change") {
		EventSampler sampler(vector<double>(1, 1.0), rng);
		while(sampler.getTime() < 1000.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
		}
		REQUIRE(sampler.getEventCount(0) > 900);
		REQUIRE(sampler.getEventCount(0) < 1100);
		
		sampler.resetEventCounts();
		
		sampler.setRate(0, 2.0);
		while(sampler.getTime() < 1500.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
		}
		REQUIRE(sampler.getEventCount(0) > 900);
		REQUIRE(sampler.getEventCount(0) < 1100);
		REQUIRE(sampler.getEventCount() == sampler.getEventCount(0));
	}
	
	SECTION("Approximately the right number of events with toggling rate.") {
		EventSampler sampler(vector<double>(1, 1.0), rng);
		bool toggle = false;
		while(sampler.getTime() < 1000.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
			sampler.setRate(0, toggle ? 2.0 : 1.0);
			toggle = !toggle;
		}
		REQUIRE(sampler.getEventCount(0) > 1230);
		REQUIRE(sampler.getEventCount(0) < 1430);
	}
	
	SECTION("Approximately the right number of events before and after several rate changes with two events.") {
		double initialRates[] = {1.0, 1.0};
		EventSampler sampler(vector<double>(begin(initialRates), end(initialRates)), rng);
		
		while(sampler.getTime() < 1000.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
		}
		REQUIRE(sampler.getEventCount(0) > 900);
		REQUIRE(sampler.getEventCount(0) < 1100);
		REQUIRE(sampler.getEventCount(1) > 900);
		REQUIRE(sampler.getEventCount(1) < 1100);
		REQUIRE(sampler.getEventCount() == sampler.getEventCount(0) + sampler.getEventCount(1));
		
		sampler.resetEventCounts();
		sampler.setRate(0, 2.0);
		
		while(sampler.getTime() < 2000.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
		}
		REQUIRE(sampler.getEventCount(0) > 1800);
		REQUIRE(sampler.getEventCount(0) < 2200);
		REQUIRE(sampler.getEventCount(1) > 900);
		REQUIRE(sampler.getEventCount(1) < 1100);
		
		sampler.resetEventCounts();
		sampler.setRate(1, 2.0);
		
		while(sampler.getTime() < 2500.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
		}
		REQUIRE(sampler.getEventCount(0) > 900);
		REQUIRE(sampler.getEventCount(0) < 1100);
		REQUIRE(sampler.getEventCount(1) > 900);
		REQUIRE(sampler.getEventCount(1) < 1100);
		
		sampler.resetEventCounts();
		sampler.setRate(0, 0.5);
		sampler.setRate(1, 0.5);
		
		while(sampler.getTime() < 4500.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
		}
		REQUIRE(sampler.getEventCount(0) > 900);
		REQUIRE(sampler.getEventCount(0) < 1100);
		REQUIRE(sampler.getEventCount(1) > 900);
		REQUIRE(sampler.getEventCount(1) < 1100);
		
		sampler.resetEventCounts();
		sampler.setRate(0, 0.5);
		sampler.setRate(1, 5.0);
		
		while(sampler.getTime() < 6500.0)
		{
			size_t eventId;
			double dt;
			sampler.sampleNextEvent(eventId, dt);
		}
		REQUIRE(sampler.getEventCount(0) > 900);
		REQUIRE(sampler.getEventCount(0) < 1100);
		REQUIRE(sampler.getEventCount(1) > 9000);
		REQUIRE(sampler.getEventCount(1) < 11000);
	}
	
	SECTION("Approximately the right number of events before and after many rate changes with many events.") {
		EventSampler sampler(vector<double>(101, 1.0), rng);
		uniform_int_distribution<size_t> ud(0, 100);
		size_t specialId = 0;
		for(int i = 0; i < 100; i++) {
			sampler.setRate(specialId, 1.0);
			specialId = ud(rng);
			sampler.setRate(specialId, 100.0);
			while(sampler.getTime() < 10.0 * (i+1)) {
				size_t eventId;
				double dt;
				sampler.sampleNextEvent(eventId, dt);
			}
			REQUIRE(sampler.getEventCount() > 1800);
			REQUIRE(sampler.getEventCount() < 2200);
			REQUIRE(sampler.getEventCount(specialId) > 900);
			REQUIRE(sampler.getEventCount(specialId) < 1100);
			sampler.resetEventCounts();
		}
	}
}
