// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

extern "C"
{
   #include <lua.h>
   #include <lauxlib.h>
   #include <lualib.h>
}

#include <iostream>
#include <map>
#include <vector>
#include <functional>
#include <cmath>
#include <limits>
#include <cassert>
#include <algorithm> 
#include <random>
#include <sstream>

using namespace std;

// #define LOG
#define TDIGEST_METATABLE "vitillo.tdigest"

struct Centroid {
public:
  Centroid(double mean, long count) : mean(mean), count(count) {}

  void add(double x, long w) {
#ifdef LOG
      cout << "Adding weight to " << this->mean << endl;
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
  friend ostream& operator<<(ostream& os, const TDigest& dt);

  TDigest(double accuracy = 0.01): m_size(0), m_accuracy(accuracy) {}

  double quantile(double q) {
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

  void add(double x, long w = 1.0) {
    vector<map<double, Centroid>::iterator> nearests = get_nearest_centroids(x);

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
      cout << "headsum " << sum << endl;
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

      // the ordering can change if the nearest point was not unique
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
      // something such as sequential ordering of data points
      // has caused a pathological expansion of our summary.
      // To fight this, we simply replay the current centroids
      // in random order.
      compress();
    }
  }

private:
  std::multimap<double, Centroid> m_centroids;
  long m_size;
  double m_accuracy;

  void add_centroid(const Centroid &centroid) {
#ifdef LOG
    cout << "Adding centroid " << centroid.mean << endl;
#endif
    m_centroids.insert({centroid.mean, centroid});
    m_size += centroid.count;
  }

  double head_sum(const map<double, Centroid>::iterator &centroid_it) {
    assert(centroid_it != m_centroids.end());

    long partial_count = 0;
    for (auto it = m_centroids.begin(); it != centroid_it; ++it) {
      partial_count += it->second.count;
    }

    return partial_count;
  }

  vector<map<double, Centroid>::iterator> get_nearest_centroids(double x) {
    vector<map<double, Centroid>::iterator> res;

    if (m_centroids.empty()) {
      return res;
    }

    auto upper_it = m_centroids.lower_bound(x);
    if (upper_it == m_centroids.begin()) {
      res.push_back(upper_it);
    } else if (upper_it == m_centroids.end()) {
      res.push_back(--upper_it);
    } else {
      auto lower_it = upper_it;
      lower_it--;
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

  void compress() {
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

};

ostream& operator<<(ostream &os, const TDigest &digest)
{
  for (const auto &pair : digest.m_centroids) {
    const Centroid &c = pair.second;
    os << "(" << c.mean << ", " << c.count << ")" << endl;
  }

  os << "Number of centroids: " << digest.m_centroids.size() << endl;
  return os;
}

static int lua_tdigest_new(lua_State *lua) {
  int n = lua_gettop(lua);
  luaL_argcheck(lua, n <= 1, 0, "incorrect number of arguments");

  double accuracy = n == 1 ? luaL_checknumber(lua, 1) : 0.01;
  size_t nbytes = sizeof(TDigest);
  TDigest **data = (TDigest **)lua_newuserdata(lua, nbytes);
  *data = new TDigest(accuracy);

  luaL_getmetatable(lua, TDIGEST_METATABLE);
  lua_setmetatable(lua, -2);

  return 1;
}

static int lua_tdigest_add(lua_State *lua) {
  int n = lua_gettop(lua);
  luaL_argcheck(lua, n == 2, 0, "incorrect number of arguments");

  TDigest **data = (TDigest **)luaL_checkudata(lua, 1, TDIGEST_METATABLE);
  double x = luaL_checknumber(lua, 2);

  (*data)->add(x);
  return 0;
}

static int lua_tdigest_quantile(lua_State *lua) {
  int n = lua_gettop(lua);
  luaL_argcheck(lua, n == 2, 0, "incorrect number of arguments");

  TDigest **data = (TDigest **)luaL_checkudata(lua, 1, TDIGEST_METATABLE);
  double q = luaL_checknumber(lua, 2);
  double x = (*data)->quantile(q);

  lua_pushnumber(lua, x);
  return 1;
}

static int lua_tdigest_tostring(lua_State *lua) {
  int n = lua_gettop(lua);
  luaL_argcheck(lua, n == 1, 0, "incorrect number of arguments");

  TDigest **data = (TDigest **)luaL_checkudata(lua, 1, TDIGEST_METATABLE);
  stringstream tmp;
  tmp << **data;
  lua_pushfstring(lua, "%s", tmp.str().c_str());
  return 1;
}

static int lua_tdigest_gc(lua_State *lua) {
  TDigest **data = (TDigest **)luaL_checkudata(lua, 1, TDIGEST_METATABLE);
  delete *data;
  return 0;
}

static const luaL_Reg tdigest_f[] = {
  {"new", lua_tdigest_new},
  {NULL, NULL}
};

static const luaL_Reg tdigest_m[] = {
  {"add", lua_tdigest_add},
  {"quantile", lua_tdigest_quantile},
  {"__tostring", lua_tdigest_tostring},
  {"__gc", lua_tdigest_gc},
  {NULL, NULL}
};

extern "C"
int luaopen_tdigest(lua_State *lua) {
  luaL_newmetatable(lua, TDIGEST_METATABLE);
  lua_pushvalue(lua, -1);  /* duplicate the metatable */
  lua_setfield(lua, -2, "__index");

  luaL_setfuncs(lua, tdigest_m, 0);
  luaL_newlib(lua, tdigest_f);
  return 1;
}

int main() {
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(1, 2);
  normal_distribution<> d1(5,1);
  normal_distribution<> d2(10,5);
  exponential_distribution<> d3(3.5);

  vector<double> sample;
  for (int i = 0; i < 10000; ++i) {
    sample.push_back(d1(gen) + d2(gen) + d3(gen));
  }

  TDigest digest;
  for (auto n : sample) {
    digest.add(n);
  }

  sort(sample.begin(), sample.end());

  for (double i = 0.1; i <= 0.99; i += 0.01) {
    double q = sample[int(i*sample.size())];
    assert(abs(q - digest.quantile(i)) < 0.1);
  }

  return 0;
}
