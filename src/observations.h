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

#ifndef OBSERVATIONS_H
#define OBSERVATIONS_H

#include <cmath>
#include <vector>
#include "global.h"

// ============================================================================
// CLASS TO CALCULATE AVG/STD + MIXTURE MODEL ON OBSERVATIONS
// ============================================================================

class Observations {

private:
        double avg;     // the average
        double M2;      // temporary variable needed for online calculation
                        // of average and variance
        size_t N;       // number of observations
        std::vector<double> observations;       // actual observations



public:
        /**
         * Default constructor
         */
        Observations() : avg(0.0), M2(0.0), N(0) {}

        /**
         * Update the average and variance but do not store the observation
         * @param x A new observation
         */
        void addObsOnline(double x) {

                N++;
                double delta = x - avg;

                avg += delta / N;
                M2 += delta*(x - avg);
        }

        const std::vector<double>& getObs() const {
                return observations;
        }

        /**
         * Store the actual observation
         * @param x A new observation
         */
        void storeObservation(double x) {
                addObsOnline(x);
                observations.push_back(x);
        }

        /**
         * Reset all parameters
         */
        void reset() {
                avg = 0.0;
                M2 = 0.0;
                N = 0;
                observations.clear();
        }

        /**
         * Get the average
         * @return The average of the observations
         */
        double getAverage() const {
                return avg;
        }

        /**
         * Get the variance
         * @return The variance
         */
        double getVariance() const {
                return M2 / N;
        }

        /**
         * Get the standard deviation
         * @return The standard deviation
         */
        double getStdev() const {
                return sqrt(getVariance());
        }

        /**
         * Get the number of observations
         * @return The number of observations
         */
        size_t getNumObs() const {
                return N;
        }
};

#endif
