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

#include "test.h"

#include "record.h"

struct Record {
  int32_t id;
  float value;
};

class RecordSerializer {
public:
  static constexpr uint32_t ENCODED_SIZE = 8;

  static void parse(Buffer& buffer, Record& record) {
    record.id = buffer.get_int32();
    record.value = buffer.get_float();
  }

  static void write(const Record& record, Buffer& buffer) {
    buffer.put_int32(record.id);
    buffer.put_float(record.value);
  }
};

static const char TEST_FILENAME[] = "RecordTest.records";

void CreateFile(size_t count) {
  FILE* file = fopen(TEST_FILENAME, "wb");
  Buffer buffer(Buffer::LITTLE, count * RecordSerializer::ENCODED_SIZE);
  for (int i = 0; i < count; ++i) {
    Record next;
    next.id = i;
    next.value = i + (float)i / count;
    RecordSerializer::write(next, buffer);
  }
  buffer.switch_mode();
  fwrite(buffer.data(), 1, buffer.remaining(), file);
  fclose(file);
}

TEST(RecordBasic) {
  Cleaner cleaner([&](){remove(TEST_FILENAME);});
  size_t count = 1000;
  CreateFile(count);
  RecordReader<Record> reader(TEST_FILENAME, &RecordSerializer::parse,
                              RecordSerializer::ENCODED_SIZE, Buffer::LITTLE, 16);
  ASSERT_EQ(0, reader.index());
  for (int i = 0; i < count; ++i) {
    reader.jump(i);
    ASSERT_EQ(i, reader.index());
    Record record;
    reader.read(record);
    ASSERT_EQ(i, record.id);
    ASSERT_EQ(i + (float)i / count, record.value, 0.000000001);
  }
  ASSERT_EQ(count, reader.index());
}

TEST(RecordJumpto) {
  Cleaner cleaner([&](){remove(TEST_FILENAME);});
  size_t count = 1000;
  CreateFile(count);
  RecordReader<Record> reader(TEST_FILENAME, &RecordSerializer::parse,
                              RecordSerializer::ENCODED_SIZE, Buffer::LITTLE, 16);
  {
    reader.jump_to([&](const Record& record)->bool{return false;});
    ASSERT_EQ(0, reader.index());
    Record record;
    reader.read(record);
    ASSERT_EQ(0, record.id);
    ASSERT_EQ(0, record.value, 0.000000001);
  }
  for (int i = 0; i < count; ++i) {
    reader.jump_to([&](const Record& record)->bool{return record.id < i;});
    ASSERT_EQ(i, reader.index());
    Record record;
    reader.read(record);
    ASSERT_EQ(i, record.id);
    ASSERT_EQ(i + (float)i / count, record.value, 0.000000001);
  }
  {
    reader.jump_to([&](const Record& record)->bool{return true;});
    ASSERT_EQ(count, reader.index());
    Record record;
    ASSERT_FALSE(reader.read(record));
  }
}
