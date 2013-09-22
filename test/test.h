/* Copyright (c) 2013 Kevin L. Stern
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#define LOG(message)\
std::cout << message << std::endl;

class TestError : public std::runtime_error {
public:
  TestError(const std::string& message) : runtime_error(message) {}
};

class TestCase {
public:
  virtual void run()=0;

  virtual const char* test_name() const=0;
};

class TestEngine {
public:
  static TestEngine& get() {
    static TestEngine instance;
    return instance;
  }

  void add_test_case(TestCase* ptr) {
    test_cases.push_back(ptr);
  }

  void run(int argc, char **argv) {
    pass_count_ = 0;
    fail_count_ = 0;
    std::unordered_set<std::string> test_cases;
    for (int i = 1; i < argc; ++i) {
      test_cases.insert(argv[i]);
    }
    if (test_cases.empty()) {
      for (const auto& next : get().test_cases) {
        LOG("Testing " << next->test_name());
        try {
          next->run();
          ++pass_count_;
        } catch (const TestError& e) {
          ++fail_count_;
          LOG("Failed: " << e.what());
        }
      }
    } else {
      for (const auto& next : get().test_cases) {
        if (test_cases.count(next->test_name()) > 0) {
          LOG("Testing " << next->test_name());
          try {
            next->run();
            ++pass_count_;
          } catch (const TestError& e) {
            ++fail_count_;
            LOG("Failed: " << e.what());
          }
        }
      }
    }
    LOG("");
    if (pass_count_ == 0 && fail_count_ == 0) {
      LOG("No tests were run.");
    } else if (fail_count_ == 0) {
      LOG("Congratulations, all tests passed!");
    } else {
      LOG(pass_count_ << "/" << (pass_count_ + fail_count_) << " tests cases passed.");
    }
  }

private:
  std::vector<TestCase*> test_cases;
  int pass_count_;
  int fail_count_;

  TestEngine() {}
  TestEngine(const TestEngine& that) {}
  TestEngine& operator=(const TestEngine& that) {return *this;}
};

#define TEST(name)\
class TEST_##name : public TestCase {\
public:\
  TEST_##name() {\
    TestEngine::get().add_test_case(this);\
  }\
  virtual ~TEST_##name() {}\
  virtual void run() {\
    run_test();\
  }\
  virtual const char* test_name() const {\
    return #name;\
  }\
private:\
  void run_test();\
};\
TEST_##name TEST_VAR_##name;\
void TEST_##name::run_test()

#define ASSERT_(oper1, oper2, op) {\
  auto r1 = oper1;\
  auto r2 = oper2;\
  if (!(r1 op r2)) {\
    throw TestError(static_cast<std::stringstream*>(&(std::stringstream() <<\
        "Expected "#oper1 << " [" << r1 << "] "#op" "#oper2 << " [" << r2 << "]"))->str());\
  }\
}

#define ASSERT_EXACT_EQ_(expected, actual)\
ASSERT_(expected, actual, ==)

#define ASSERT_FUZZY_EQ_(expected, actual, tolerance)\
ASSERT_(abs(expected - actual), tolerance, <)

#define GET_ASSERT_MACRO_(p1, p2, p3, name, ...) name

#define ASSERT_EQ(...)\
GET_ASSERT_MACRO_(__VA_ARGS__, ASSERT_FUZZY_EQ_, ASSERT_EXACT_EQ_)(__VA_ARGS__)

#define ASSERT_ARRAY_EQ(expected, actual, size) {\
  bool pass = true;\
  for (int i = 0; i < size; ++i) {\
    if (expected[i] != actual[i]) {\
      pass = false;\
    }\
  }\
  if (!pass) {\
    std::stringstream s;\
    s << "[";\
    for (int i = 0; i < size; ++i) {\
      s << expected[i];\
      if (i < size - 1) {\
        s << ",";\
      }\
    }\
    s << "]";\
    s << " != ";\
    s << "[";\
    for (int i = 0; i < size; ++i) {\
      s << actual[i];\
      if (i < size - 1) {\
        s << ",";\
      }\
    }\
    s << "]";\
    throw TestError(s.str());\
  }\
}

#define ASSERT_NEQ(expected, actual)\
ASSERT_(expected, actual, !=)

#define ASSERT_LT(smaller, larger)\
ASSERT_(smaller, larger, <)

#define ASSERT_LTE(not_larger, not_smaller)\
ASSERT_(not_larger, not_smaller, <=)

#define ASSERT_GT(larger, smaller)\
ASSERT_(larger, smaller, >)

#define ASSERT_GTE(not_smaller, not_larger)\
ASSERT_(not_smaller, not_larger, >=)

#define ASSERT_TRUE(statement)\
if (!(statement)) {\
  throw TestError("Expected "#statement" [false] true");\
}

#define ASSERT_FALSE(statement)\
if (statement) {\
  throw TestError("Expected "#statement" [true] false");\
}

#define ASSERT_NULL(statement)\
if (statement != nullptr) {\
  throw TestError("Expected "#statement" null");\
}

// Utility class for cleaning a resource when an instance falls out of scope.
class Cleaner {
public:
  Cleaner(std::function<void()> cleanup) : cleanup_(cleanup) {}

  ~Cleaner() {
    cleanup_();
  }

private:
  std::function<void()> cleanup_;
};
