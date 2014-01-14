/* Copyright (c) 2014 Kevin L. Stern
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

#include "test.h"

#include "linear_hashing_table.h"

namespace {
static const char TEST_PRIMARY_FILENAME[] = "LinearHashingTableTest.primary_records";
static const char TEST_OVERFLOW_FILENAME[] = "LinearHashingTableTest.overflow_records";

class RecordSerializer {
public:
  static constexpr uint32_t ENCODED_SIZE = 8;

  static void parse(Buffer& buffer, int32_t& key, float& value) {
    key = buffer.get_int32();
    value = buffer.get_float();
  }

  static void write(const int32_t key, const float value, Buffer& buffer) {
    buffer.put_int32(key);
    buffer.put_float(value);
  }
};

}

TEST(LinearHashingTableInMemory) {
  const int iters = 1000;
  LinearHashingTable<int32_t, float>::TransientBucketManager manager(10);
  LinearHashingTable<int32_t, float> table(&manager);
  for (int i = 0; i < iters; ++i) {
    ASSERT_FALSE(table.contains(i));
    table.insert(i, static_cast<float>(i) / iters);
    for (int j = 0; j <= i; ++j) {
      ASSERT_TRUE(table.contains(j));
    }
  }
  ASSERT_EQ(iters, table.size());
  float value;
  for (int i = iters - 1; i >= 0; --i) {
    ASSERT_TRUE(table.remove(i, value));
    ASSERT_EQ(static_cast<float>(i) / iters, value);
    ASSERT_FALSE(table.contains(i));
    for (int j = 0; j < i; ++j) {
      ASSERT_TRUE(table.contains(j));
    }
  }
}

#define INSTANTIATE_TABLE \
LinearHashingTable<int32_t, float>::FileBucketManager manager( \
    10 /* bucket capacity */, TEST_PRIMARY_FILENAME, TEST_OVERFLOW_FILENAME, \
    &RecordSerializer::write, &RecordSerializer::parse, RecordSerializer::ENCODED_SIZE, \
    Buffer::LITTLE, 10 /* bucket cache size */); \
LinearHashingTable<int32_t, float> table(&manager);

TEST(LinearHashingTableFilesystem) {
  Cleaner cleaner1([&](){remove(TEST_PRIMARY_FILENAME);});
  Cleaner cleaner2([&](){remove(TEST_OVERFLOW_FILENAME);});
  const int iters = 100;

  const int lapcount = 10;
  const int iters_per_lap = iters / lapcount;
  for (int lap = 0; lap < lapcount; ++lap) {
    INSTANTIATE_TABLE;
    size_t initial_table_size = table.size();
    for (int i = 0; i < (lap + 1) * iters_per_lap; ++i) {
      if (i >= lap * iters_per_lap) {
        ASSERT_TRUE(table.insert(i, static_cast<float>(i) / iters));
        ASSERT_EQ(i + 1, table.size());
      } else {
        ASSERT_FALSE(table.insert(i, static_cast<float>(i) / iters));
        ASSERT_EQ(initial_table_size, table.size());
      }
      for (int j = 0; j < i; ++j) {
        ASSERT_TRUE(table.contains(j));
      }
    }
  }
  for (int lap = 0; lap < lapcount; ++lap) {
    INSTANTIATE_TABLE;
    size_t initial_table_size = table.size();
    float value;
    for (int i = 0; i < (lap + 1) * iters_per_lap; ++i) {
      if (i >= lap * iters_per_lap) {
        ASSERT_TRUE(table.find(i, value));
        ASSERT_EQ(static_cast<float>(i) / iters, value);
        ASSERT_TRUE(table.remove(i, value));
        ASSERT_EQ(static_cast<float>(i) / iters, value);
        ASSERT_EQ(iters - (i + 1), table.size());
      } else {
        ASSERT_FALSE(table.remove(i, value));
        ASSERT_EQ(initial_table_size, table.size());
      }
      for (int j = 0; j < i; ++j) {
        ASSERT_FALSE(table.contains(j));
      }
    }
  }
}

#undef INSTANTIATE_TABLE
