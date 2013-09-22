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

#include "buffer.h"

TEST(BufferLittleEndian) {
  {
    uint8_t bytes[4];
    for (int i = 0; i < 4; ++i) {
      bytes[i] = i;
    }
    int32_t val = *reinterpret_cast<int32_t*>(bytes);
    Buffer buffer(Buffer::LITTLE, 32);
    buffer.put_int32(val);
    buffer.switch_mode();
    buffer.get_bytes(bytes, 4);
    for (int i = 0; i < 4; ++i) {
      ASSERT_EQ(i, bytes[i]);
    }
    buffer.reset();
    ASSERT_EQ(val, buffer.get_int32());
  }
  {
    uint8_t bytes[4];
    for (int i = 0; i < 4; ++i) {
      bytes[i] = i;
    }
    uint32_t val = *reinterpret_cast<uint32_t*>(bytes);
    Buffer buffer(Buffer::LITTLE, 32);
    buffer.put_uint32(val);
    buffer.switch_mode();
    buffer.get_bytes(bytes, 4);
    for (int i = 0; i < 4; ++i) {
      ASSERT_EQ(i, bytes[i]);
    }
    buffer.reset();
    ASSERT_EQ(val, buffer.get_uint32());
  }
  {
    uint8_t bytes[8];
    for (int i = 0; i < 8; ++i) {
      bytes[i] = i;
    }
    int64_t val = *reinterpret_cast<int64_t*>(bytes);
    Buffer buffer(Buffer::LITTLE, 64);
    buffer.put_int64(val);
    buffer.switch_mode();
    buffer.get_bytes(bytes, 8);
    for (int i = 0; i < 8; ++i) {
      ASSERT_EQ(i, bytes[i]);
    }
    buffer.reset();
    ASSERT_EQ(val, buffer.get_int64());
  }
  {
    uint8_t bytes[8];
    for (int i = 0; i < 8; ++i) {
      bytes[i] = i;
    }
    uint64_t val = *reinterpret_cast<uint64_t*>(bytes);
    Buffer buffer(Buffer::LITTLE, 64);
    buffer.put_uint64(val);
    buffer.switch_mode();
    buffer.get_bytes(bytes, 8);
    for (int i = 0; i < 8; ++i) {
      ASSERT_EQ(i, bytes[i]);
    }
    buffer.reset();
    ASSERT_EQ(val, buffer.get_uint64());
  }
}

TEST(BufferBigEndian) {
  {
    uint8_t bytes[4];
    for (int i = 0; i < 4; ++i) {
      bytes[i] = i;
    }
    int32_t val = *reinterpret_cast<int32_t*>(bytes);
    Buffer buffer(Buffer::BIG, 32);
    buffer.put_int32(val);
    buffer.switch_mode();
    buffer.get_bytes(bytes, 4);
    for (int i = 0; i < 4; ++i) {
      ASSERT_EQ(i, bytes[3-i]);
    }
    buffer.reset();
    ASSERT_EQ(val, buffer.get_int32());
  }
  {
    uint8_t bytes[4];
    for (int i = 0; i < 4; ++i) {
      bytes[i] = i;
    }
    uint32_t val = *reinterpret_cast<uint32_t*>(bytes);
    Buffer buffer(Buffer::BIG, 32);
    buffer.put_uint32(val);
    buffer.switch_mode();
    buffer.get_bytes(bytes, 4);
    for (int i = 0; i < 4; ++i) {
      ASSERT_EQ(i, bytes[3-i]);
    }
    buffer.reset();
    ASSERT_EQ(val, buffer.get_uint32());
  }
  {
    uint8_t bytes[8];
    for (int i = 0; i < 8; ++i) {
      bytes[i] = i;
    }
    int64_t val = *reinterpret_cast<int64_t*>(bytes);
    Buffer buffer(Buffer::BIG, 64);
    buffer.put_int64(val);
    buffer.switch_mode();
    buffer.get_bytes(bytes, 8);
    for (int i = 0; i < 8; ++i) {
      ASSERT_EQ(i, bytes[7-i]);
    }
    buffer.reset();
    ASSERT_EQ(val, buffer.get_int64());
  }
  {
    uint8_t bytes[8];
    for (int i = 0; i < 8; ++i) {
      bytes[i] = i;
    }
    uint64_t val = *reinterpret_cast<uint64_t*>(bytes);
    Buffer buffer(Buffer::BIG, 64);
    buffer.put_uint64(val);
    buffer.switch_mode();
    buffer.get_bytes(bytes, 8);
    for (int i = 0; i < 8; ++i) {
      ASSERT_EQ(i, bytes[7-i]);
    }
    buffer.reset();
    ASSERT_EQ(val, buffer.get_uint64());
  }
}

TEST(BufferDrainTo) {
  const char* phrase = "Hello World!";
  size_t phrase_length = strlen(phrase);
  Buffer buffer(phrase_length);
  for (size_t i = 0; i < phrase_length; ++i) {
    buffer.put_byte(phrase[i]);
  }
  buffer.switch_mode();
  std::stringstream stream;
  buffer.drain_to(stream);
  ASSERT_EQ(phrase, stream.str());
}
