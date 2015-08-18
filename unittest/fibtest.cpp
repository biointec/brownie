#include <gtest/gtest.h>
#include <iostream>

using namespace std;

#include "fiblib/fibheap.h"

TEST(FibHeap, FibHeapTest)
{
        FibHeap<int> fh;
        FibHeapEl *el[10];

        for (int i = 0; i < 10; i++)
                el[i] = fh.insert(0, i);

        fh.replaceKey(el[1], -1);
        fh.replaceKey(el[6], -1);
        fh.replaceKey(el[4], -1);
        fh.replaceKey(el[2], -1);
        fh.replaceKey(el[8], -1);

        ASSERT_EQ(fh.getMinKey(), -1);
        ASSERT_EQ(fh.extractMinKey(), 8);

        fh.replaceKey(el[7], -33);
        fh.replaceKey(el[4], -36);
        fh.replaceKey(el[3], -1);
        fh.replaceKey(el[9], -81.2);

        ASSERT_FLOAT_EQ(fh.getMinKey(), -81.2);
        ASSERT_EQ(fh.extractMinKey(), 9);

        fh.replaceKey(el[6], -68);
        fh.replaceKey(el[2], -69);

        ASSERT_EQ(fh.getMinKey(), -69);
        ASSERT_EQ(fh.extractMinKey(), 2);

        fh.replaceKey(el[1], -52);
        fh.replaceKey(el[3], -2);
        fh.replaceKey(el[4], -120);
        fh.replaceKey(el[5], -48);

        ASSERT_EQ(fh.getMinKey(), -120);
        ASSERT_EQ(fh.extractMinKey(), 4);

        fh.replaceKey(el[3], -3);
        fh.replaceKey(el[5], -63);

        ASSERT_EQ(fh.getMinKey(), -68);
        ASSERT_EQ(fh.extractMinKey(), 6);

        fh.replaceKey(el[5], -110);
        fh.replaceKey(el[7], -115);

        ASSERT_EQ(fh.getMinKey(), -115);
        ASSERT_EQ(fh.extractMinKey(), 7);

        fh.replaceKey(el[5], -188);

        ASSERT_EQ(fh.getMinKey(), -188);
        ASSERT_EQ(fh.extractMinKey(), 5);
}
