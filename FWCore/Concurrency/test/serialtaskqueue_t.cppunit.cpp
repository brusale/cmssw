//
//  SerialTaskQueue_test.cpp
//  DispatchProcessingDemo
//
//  Created by Chris Jones on 9/27/11.
//

#include <iostream>

#include <cppunit/extensions/HelperMacros.h>
#include <chrono>
#include <memory>
#include <atomic>
#include "oneapi/tbb/task_arena.h"
#include "FWCore/Concurrency/interface/WaitingTask.h"
#include "FWCore/Concurrency/interface/SerialTaskQueue.h"
#include "FWCore/Concurrency/interface/FunctorTask.h"

class SerialTaskQueue_test : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(SerialTaskQueue_test);
  CPPUNIT_TEST(testPush);
  CPPUNIT_TEST(testPause);
  CPPUNIT_TEST(stressTest);
  CPPUNIT_TEST_SUITE_END();

public:
  void testPush();
  void testPause();
  void stressTest();
  void setUp() {}
  void tearDown() {}
};

CPPUNIT_TEST_SUITE_REGISTRATION(SerialTaskQueue_test);
using namespace std::chrono_literals;

void SerialTaskQueue_test::testPush() {
  std::atomic<unsigned int> count{0};

  edm::SerialTaskQueue queue;
  {
    std::atomic<unsigned int> waitingTasks{3};
    oneapi::tbb::task_group group;

    queue.push(group, [&count, &waitingTasks] {
      CPPUNIT_ASSERT(count++ == 0);
      std::this_thread::sleep_for(10us);
      --waitingTasks;
    });

    queue.push(group, [&count, &waitingTasks] {
      CPPUNIT_ASSERT(count++ == 1);
      std::this_thread::sleep_for(10us);
      --waitingTasks;
    });

    queue.push(group, [&count, &waitingTasks] {
      CPPUNIT_ASSERT(count++ == 2);
      std::this_thread::sleep_for(10us);
      --waitingTasks;
    });

    do {
      group.wait();
    } while (0 != waitingTasks.load());
    CPPUNIT_ASSERT(count == 3);
  }
}

void SerialTaskQueue_test::testPause() {
  std::atomic<unsigned int> count{0};

  edm::SerialTaskQueue queue;
  {
    queue.pause();
    {
      std::atomic<unsigned int> waitingTasks{1};
      oneapi::tbb::task_group group;

      queue.push(group, [&count, &waitingTasks] {
        CPPUNIT_ASSERT(count++ == 0);
        --waitingTasks;
      });
      std::this_thread::sleep_for(1000us);
      CPPUNIT_ASSERT(0 == count);
      queue.resume();
      do {
        group.wait();
      } while (0 != waitingTasks.load());
      CPPUNIT_ASSERT(count == 1);
    }

    {
      std::atomic<unsigned int> waitingTasks{3};
      oneapi::tbb::task_group group;

      queue.push(group, [&count, &queue, &waitingTasks] {
        queue.pause();
        CPPUNIT_ASSERT(count++ == 1);
        --waitingTasks;
      });
      queue.push(group, [&count, &waitingTasks] {
        CPPUNIT_ASSERT(count++ == 2);
        --waitingTasks;
      });
      queue.push(group, [&count, &waitingTasks] {
        CPPUNIT_ASSERT(count++ == 3);
        --waitingTasks;
      });
      std::this_thread::sleep_for(100us);
      //can't do == since the queue may not have processed the first task yet
      CPPUNIT_ASSERT(2 >= count);
      queue.resume();
      do {
        group.wait();
      } while (0 != waitingTasks.load());
      CPPUNIT_ASSERT(count == 4);
    }
  }
}

void SerialTaskQueue_test::stressTest() {
  //note group needs to live longer than queue
  oneapi::tbb::task_group group;
  edm::SerialTaskQueue queue;

  unsigned int index = 100;
  const unsigned int nTasks = 1000;
  while (0 != --index) {
    std::atomic<unsigned int> waitingTasks{2};
    std::atomic<unsigned int> count{0};

    std::atomic<bool> waitToStart{true};
    {
      group.run([&queue, &waitToStart, &waitingTasks, &count, &group] {
        while (waitToStart.load()) {
        }
        for (unsigned int i = 0; i < nTasks; ++i) {
          ++waitingTasks;
          queue.push(group, [&count, &waitingTasks] {
            ++count;
            --waitingTasks;
          });
        }

        --waitingTasks;
      });

      group.run([&queue, &waitToStart, &waitingTasks, &count, &group] {
        waitToStart = false;
        for (unsigned int i = 0; i < nTasks; ++i) {
          ++waitingTasks;
          queue.push(group, [&count, &waitingTasks] {
            ++count;
            --waitingTasks;
          });
        }
        --waitingTasks;
      });
    }
    do {
      group.wait();
    } while (0 != waitingTasks.load());

    CPPUNIT_ASSERT(2 * nTasks == count);
  }
}

#include <Utilities/Testing/interface/CppUnit_testdriver.icpp>
