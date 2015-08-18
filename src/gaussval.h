/***************************************************************************
 *   Copyright (C) 2014, 2015 Jan Fostier (jan.fostier@intec.ugent.be)     *
 *   Copyright (C) 2014, 2015 Mahdi Heydari (mahdi.heydari@intec.ugent.be) *
 *   This file is part of Brownie                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef GAUSSVAL_H
#define GAUSSVAL_H

#define DSMALL 1E-13
#define DSIGMA 3

#include <cmath>

// ============================================================================
// GAUSSIAN VALUE (MU, SIGMA) CLASS
// ============================================================================

class GaussVal
{
private:
        double avg;     // average mu
        double std;     // standard deviation

public:
        /**
         * Default constructor
         */
        GaussVal() : avg(0.0), std(0.0) {}

        /**
         * Constructor
         * @param avg Average value
         * @param std Standard deviation
         */
        GaussVal(double avg_, double std_) : avg(avg_), std(std_) {}

        /**
         * Set the average value
         * @param target Target value
         */
        void setAverage(double target) {
                avg = target;
        }

        /**
         * Get the average value
         * @return The average value
         */
        double getAverage() const {
                return avg;
        }

        /**
         * Set the standard deviation
         * @param target Target value
         */
        void setStd(double target) {
                std = target;
        }

        /**
         * Get the standard deviation
         * @return The standard deviation
         */
        double getStd() const {
                return std;
        }

        /**
         * Check whether this gauss value is significantly smaller than the rhs
         * @param rhs Right hand side
         * @return True of false
         */
        bool isSignSmaller(const GaussVal& rhs) const {
                if (avg + DSIGMA * std > rhs.getAverage())
                        return false;
                if (avg > rhs.getAverage() - DSIGMA * rhs.getStd())
                        return false;
                return true;
        }

        /**
         * Check whether this gauss value is significantly bigger than the rhs
         * @param rhs Right hand side
         * @return True of false
         */
        bool isSignBigger(const GaussVal& rhs) const {
                if (avg - DSIGMA * std < rhs.getAverage())
                        return false;
                if (avg < rhs.getAverage() + DSIGMA * rhs.getStd())
                        return false;
                return true;
        }

         /**
         * Check whether this gauss value is significantly smaller than the rhs
         * @param rhs Right hand side
         * @return True of false
         */
        bool isSignSmaller(double rhs) const {
                if (avg + DSIGMA * std > rhs)
                        return false;
                return true;
        }

        /**
         * Check whether this gauss value is significantly bigger than the rhs
         * @param rhs Right hand side
         * @return True of false
         */
        bool isSignBigger(double rhs) const {
                if (avg - DSIGMA * std < rhs)
                        return false;
                return true;
        }

        /**
         * Check whether two gauss intervals overlaps
         * @param leftA left gauss value for interval A
         * @param rightA right gauss value for interval A
         * @param leftB left gauss value for interval B
         * @param rightB right gauss value for interval B
         */
        static bool overlaps(const GaussVal& leftA, const GaussVal& rightA,
                             const GaussVal& leftB, const GaussVal& rightB)
        {
                if (rightA.isSignSmaller(leftB))
                        return false;
                if (rightB.isSignSmaller(leftA))
                        return false;
                return true;
        }

        /**
         * Operator < overloading
         * @param rhs A const-reference to the right-hand-side
         * @return True if average is smaller than rhs.average
         */
        bool operator<(const GaussVal& rhs) const {
                return avg < rhs.avg;
        }

        bool operator==(const GaussVal& rhs) const {
                return ((std::abs(avg - rhs.avg) < DSMALL) &&
                        (std::abs(std - rhs.std) < DSMALL));
        }

        bool operator!=(const GaussVal& rhs) const {
                return !(*this == rhs);
        }

        /**
         * Operator + overloading
         * @param rhs A const-reference to the right-hand side
         * @return A Gaussval containing the sum of the two objects
         */
        GaussVal operator+(const GaussVal& rhs) const {
                return GaussVal(avg + rhs.avg, sqrt(std*std + rhs.std*rhs.std));
        }

        /**
         * Operator - overloading
         * @param rhs A const-reference to the right-hand side
         * @return A Gaussval containing the subtraction of the two objects
         */
        GaussVal operator-(const GaussVal& rhs) const {
                return GaussVal(avg - rhs.avg, sqrt(std*std + rhs.std*rhs.std));
        }

        /**
         * Operator + overloading
         * @param rhs A const-reference to the right-hand side
         * @return A Gaussval containing the sum of the two objects
         */
        GaussVal operator+(double rhs) const {
                return GaussVal(avg + rhs, std);
        }

        /**
         * Operator + overloading
         * @param rhs A const-reference to the right-hand side
         * @return A Gaussval containing the sum of the two objects
         */
        GaussVal operator-(double rhs) const {
                return GaussVal(avg - rhs, std);
        }

        /**
         * Operator * overloading
         * @param rhs A const-reference to the right-hand side
         * @return A Gaussval containing the product of the two objects
         */
        GaussVal operator*(double rhs) const {
                return GaussVal(avg * rhs, std * std::abs(rhs));
        }

        friend std::ostream &operator<<(std::ostream &out, const GaussVal &gv) {
               // out << gv.getAverage() << " +/- " << gv.getStd();
                return out;
        }
};

class SortGVByStd {

public:
        bool operator()(const GaussVal& lhs, const GaussVal& rhs) {
                return lhs.getStd() < rhs.getStd();
        }
};

#endif
