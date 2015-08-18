#include "readfile/fastafile.h"
#include "readfile/fastqfile.h"
#include "readfile/rawfile.h"
#include "readfile/samfile.h"

#include <gtest/gtest.h>
#include <cstdio>
#include <string>

using namespace std;

// ============================================================================
// SAM FORMAT
// ============================================================================

TEST(readFile, SAMGZTest)
{
        ReadFile *file = new SamFile(true);
        file->open("test.sam.gz");

        string target("TAATCCCCGCCAAATTCGTGACCTGTCATTCGTCG");
        string read, description;
        file->getNextRead(read, description);

        EXPECT_EQ(read, target);

        target = "CGACAGGGATAGTGTAGCTGACCGTTGTGACTGGC";

        int numReads = 1;
        while (file->getNextRead(read, description)) {
                numReads++;

                if (numReads == 10)
                        EXPECT_EQ(read, target);
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}

TEST(readFile, SAMTest)
{
        ReadFile *file = new SamFile(false);
        file->open("test.sam");

        string target("TAATCCCCGCCAAATTCGTGACCTGTCATTCGTCG");
        string read, description;
        file->getNextRead(read, description);

        EXPECT_EQ(read, target);

        target = "CGACAGGGATAGTGTAGCTGACCGTTGTGACTGGC";

        int numReads = 1;
        while (file->getNextRead(read, description)) {
                numReads++;

                if (numReads == 10)
                        EXPECT_EQ(read, target);
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}

// ============================================================================
// FASTA FORMAT
// ============================================================================

TEST(readFile, FastAGZTest)
{
        ReadFile *file = new FastAFile(true);
        file->open("test.fa.gz");

        string target("GAACGGTCCGGCCGCATCCATTTCTTCCCTGTAGCGAATCGCGAAAATCGTCCGGA"
                      "GTCTTAGTGTCTAAAGGTGGTTCACACGGAGATATGAGCGCGCC");
        string read, description;
        file->getNextRead(read, description);

        EXPECT_EQ(read, target);

        target = "GGTGAACTGACTCCGAGTGCCATACTGTTCTTCGACCACGGTCATGGACCCAGTCGTATT"
                 "TAGAAATGCCCCAGGGAAGAGGTATAGTAGATAGACGGCT";

        int numReads = 1;
        while (file->getNextRead(read, description)) {
                numReads++;

                if (numReads == 10)
                        EXPECT_EQ(read, target);
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}

TEST(readFile, FastATest)
{
        ReadFile *file = new FastAFile(false);
        file->open("test.fa");

        string target("GAACGGTCCGGCCGCATCCATTTCTTCCCTGTAGCGAATCGCGAAAATCGTCCGGA"
                      "GTCTTAGTGTCTAAAGGTGGTTCACACGGAGATATGAGCGCGCC");
        string read, description;
        file->getNextRead(read, description);

        EXPECT_EQ(read, target);

        target = "GGTGAACTGACTCCGAGTGCCATACTGTTCTTCGACCACGGTCATGGACCCAGTCGTATT"
                 "TAGAAATGCCCCAGGGAAGAGGTATAGTAGATAGACGGCT";

        int numReads = 1;
        while (file->getNextRead(read, description)) {
                numReads++;

                if (numReads == 10)
                        EXPECT_EQ(read, target);
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}

// ============================================================================
// FASTQ FORMAT
// ============================================================================

TEST(readFile, FastQTest)
{
        ReadFile *file = new FastQFile(true);
        file->open("test.fastq.gz");

        string target("ATATAGATGTACATAAATTAGTTGAAGTATATGAACG");
        string read, description;
        file->getNextRead(read, description);

        EXPECT_EQ(read, target);

        target = "TAGGAAAGCGAAGCCATTCAATACGAAGTATTGTATA";

        int numReads = 1;
        while (file->getNextRead(read, description)) {
                numReads++;

                if (numReads == 9)
                        EXPECT_EQ(read, target);
        }

        EXPECT_EQ(numReads, 9);

        file->close();
        delete file;
}


TEST(readFile, FastQGZTest)
{
        ReadFile *file = new FastQFile(false);
        file->open("test.fastq");

        string target("ATATAGATGTACATAAATTAGTTGAAGTATATGAACG");
        string read, description;
        file->getNextRead(read, description);

        EXPECT_EQ(read, target);

        target = "TAGGAAAGCGAAGCCATTCAATACGAAGTATTGTATA";

        int numReads = 1;
        while (file->getNextRead(read, description)) {
                numReads++;

                if (numReads == 9)
                        EXPECT_EQ(read, target);
        }

        EXPECT_EQ(numReads, 9);

        file->close();
        delete file;
}

// ============================================================================
// RAW FORMAT
// ============================================================================

TEST(readFile, RawGZTest)
{
        ReadFile *file = new RawFile(true);
        file->open("test.raw.gz");

        string target("GAACGGTCCGGCCGCATCCATTTCTTCCCTGTAGCGAATCGCGAAAATCGTCCGGA"
                      "GTCTTAGTGTCTAAAGGTGGTTCACACGGAGATATGAGCGCGCC");
        string read, description;
        file->getNextRead(read, description);

        EXPECT_EQ(read, target);

        target = "GGTGAACTGACTCCGAGTGCCATACTGTTCTTCGACCACGGTCATGGACCCAGTCGTATT"
                 "TAGAAATGCCCCAGGGAAGAGGTATAGTAGATAGACGGCT";

        int numReads = 1;
        while (file->getNextRead(read, description)) {
                numReads++;

                if (numReads == 10)
                        EXPECT_EQ(read, target);
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}

TEST(readFile, RawTest)
{
        ReadFile *file = new RawFile(false);
        file->open("test.raw");

        string target("GAACGGTCCGGCCGCATCCATTTCTTCCCTGTAGCGAATCGCGAAAATCGTCCGGA"
                      "GTCTTAGTGTCTAAAGGTGGTTCACACGGAGATATGAGCGCGCC");

        string read, description;
        file->getNextRead(read, description);

        EXPECT_EQ(read, target);

        target = "GGTGAACTGACTCCGAGTGCCATACTGTTCTTCGACCACGGTCATGGACCCAGTCGTATT"
                 "TAGAAATGCCCCAGGGAAGAGGTATAGTAGATAGACGGCT";

        int numReads = 1;
        while (file->getNextRead(read, description)) {
                numReads++;

                if (numReads == 10)
                        EXPECT_EQ(read, target);
        }

        EXPECT_EQ(numReads, 10);

        file->close();
        delete file;
}
