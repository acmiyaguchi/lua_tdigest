// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <random>

#include "TDigest.h"

using namespace std;

void TDigest::add(double x, long w) {
  auto nearests = get_nearest_centroids(x);

#ifdef LOG
  cout << "Nearest neighbours for " << x << endl;
  for (auto c : nearests)
    cout << c->second.mean << " - " << c->second.count << endl;
#endif

  if (nearests.empty()) {
    add_centroid(Centroid(x, w));
    return;
  }

  long sum = head_sum(nearests[0]);
  random_shuffle(nearests.begin(), nearests.end());

#ifdef LOG
  cout << "headsum: " << sum << endl;
#endif

  for (size_t i = 0; i < nearests.size() && w > 0; ++i) {
    Centroid centroid  = nearests[i]->second;
    double q = (sum + centroid.count*0.5)/m_size;
    double limit = 4*m_size*q*(1-q)*m_accuracy;

#ifdef LOG
    cout << "limit: " << limit << endl;
#endif

    if (centroid.count > limit) {
      continue;
    }

    // The ordering can change if the nearest point isn't unique.
    long dw = min(long(limit - centroid.count), w);
    m_size -= centroid.count; //TODO: single op in erase
    m_centroids.erase(nearests[i]);
    centroid.add(x, dw);
    add_centroid(centroid);

    w -= dw;
    sum += centroid.count;
  }

  if (w > 0) {
    add_centroid(Centroid(x, w));
  }

  if (m_centroids.size() > 10 * 1/m_accuracy) {
    // Something such as sequential ordering of data points has caused a
    // pathological expansion of our summary. To fight this, we simply
    // replay the current centroids in random order.
    compress();
  }
}

double TDigest::quantile(double q) {
  assert(q >= 0 && q <= 1);

  if (m_centroids.empty()) {
    return nan("");
  }

  if (q == 0) {
    return 0;
  }

  if(m_centroids.size() == 1) {
    return m_centroids.begin()->second.mean;
  }

  q *= m_size;
  long sum = 0;

  for (auto it = m_centroids.begin(); it != m_centroids.end(); ++it) {
    const Centroid &c = it->second;

    if (q < sum + c.count) {
      double delta = 0;

      if (it == m_centroids.begin()) {
        delta = (++it)->second.mean - c.mean;
      } else if (it == m_centroids.end()) {
        delta = c.mean - (--it)->second.mean;
      } else {
        auto after = it; after++;
        delta = 0.5*(after->second.mean - (--it)->second.mean);
      }

      return c.mean + (double(q - sum)/c.count - 0.5)*delta;
    }

    sum += c.count;
  }

  return m_centroids.rbegin()->second.mean;
}

void TDigest::add_centroid(const Centroid &centroid) {
#ifdef LOG
  cout << "Adding centroid " << centroid.mean << endl;
#endif
  m_centroids.insert({centroid.mean, centroid});
  m_size += centroid.count;
}

// Note, this method could be made faster using a smarter data structure.
double TDigest::head_sum(const map<double, Centroid>::iterator &centroid_it) {
  assert(centroid_it != m_centroids.end());

  long partial_count = 0;
  for (auto it = m_centroids.begin(); it != centroid_it; ++it) {
    partial_count += it->second.count;
  }

  return partial_count;
}

vector<TDigest::CentroidMapType::iterator> TDigest::get_nearest_centroids(double x) {
  vector<CentroidMapType::iterator> res;

  if (m_centroids.empty()) {
    return res;
  }

  auto upper_it = m_centroids.lower_bound(x);

  if (upper_it == m_centroids.begin()) {
    res.push_back(upper_it);
  } else if (upper_it == m_centroids.end()) {
    res.push_back(--upper_it);
  } else {
    auto lower_it = upper_it; lower_it--;
    Centroid upper = upper_it->second;
    Centroid lower = lower_it->second;
    double distance_diff = lower.distance(x) - upper.distance(x);

    if (distance_diff <= 0) {
      for (; lower_it->first == lower.mean && lower_it != m_centroids.begin(); --lower_it) {
        res.push_back(lower_it);
      }
      reverse(res.begin(), res.end());
    }
    if (distance_diff >= 0) {
      for(; upper_it->first == upper.mean && upper_it != m_centroids.end(); ++upper_it) {
        res.push_back(upper_it);
      }
    }
  }

  return res;
}

void TDigest::compress() {
  TDigest compressed(m_accuracy);
  vector<Centroid> tmp;

#ifdef LOG
  cout << "Compressing..." << endl;
#endif

  for (const auto &it : m_centroids) {
    tmp.push_back(it.second);
  }

  random_shuffle(tmp.begin(), tmp.end());
  for (const auto &c : tmp) {
    compressed.add(c.mean, c.count);
  }

  m_centroids = move(compressed.m_centroids);
}

ostream& operator<<(ostream &os, const TDigest &digest)
{
  for (const auto &pair : digest.m_centroids) {
    const Centroid &c = pair.second;
    os << "(" << c.mean << ", " << c.count << ")" << endl;
  }

  os << "Number of centroids: " << digest.m_centroids.size() << endl;
  return os;
}
