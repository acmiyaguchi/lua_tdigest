// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef TDIGEST_H
#define TDIGEST_H

#include <iostream>
#include <map>
#include <vector>

struct Centroid {
public:
  Centroid(double mean, long count) : mean(mean), count(count) {}

  void add(double x, long w) {
#ifdef LOG
    std::cout << "Adding weight to " << this->mean << endl;
#endif
    this->count += w;
    this->mean += w * (x - this->mean)/count;
  }

  double distance(const Centroid &other) {
    return abs(this->mean - other.mean);
  }

  double distance(double mean) {
    return abs(this->mean - mean);
  }

  double mean;
  long count;
};

class TDigest {
public:
  friend std::ostream& operator<<(std::ostream& os, const TDigest& dt);

  TDigest(double accuracy = 0.01): m_size(0), m_accuracy(accuracy) {}

  void add(double x, long w = 1.0);
  double quantile(double q);

 private:
  typedef typename std::multimap<double, Centroid> CentroidMapType;

  CentroidMapType m_centroids;
  long m_size;
  double m_accuracy;

  void add_centroid(const Centroid &centroid);
  double head_sum(const CentroidMapType::iterator &centroid_it);
  std::vector<CentroidMapType::iterator> get_nearest_centroids(double x);
  void compress();
};

#endif
