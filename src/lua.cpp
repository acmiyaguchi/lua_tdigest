// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sstream>

#include "TDigest.h"

using namespace std;

extern "C"
{
   #include <lua.h>
   #include <lauxlib.h>
   #include <lualib.h>
}

const char *tdigest_metatable = "vitillo.tdigest";

static int lua_tdigest_new(lua_State *lua) {
  int n = lua_gettop(lua);
  luaL_argcheck(lua, n <= 1, 0, "incorrect number of arguments");

  double accuracy = n == 1 ? luaL_checknumber(lua, 1) : 0.01;
  size_t nbytes = sizeof(TDigest*);
  TDigest **data = (TDigest **)lua_newuserdata(lua, nbytes);
  *data = new TDigest(accuracy);

  luaL_getmetatable(lua, tdigest_metatable);
  lua_setmetatable(lua, -2);

  return 1;
}

static int lua_tdigest_add(lua_State *lua) {
  int n = lua_gettop(lua);
  luaL_argcheck(lua, n == 2, 0, "incorrect number of arguments");

  TDigest **data = (TDigest **)luaL_checkudata(lua, 1, tdigest_metatable);
  double x = luaL_checknumber(lua, 2);

  (*data)->add(x);
  return 0;
}

static int lua_tdigest_quantile(lua_State *lua) {
  int n = lua_gettop(lua);
  luaL_argcheck(lua, n == 2, 0, "incorrect number of arguments");

  TDigest **data = (TDigest **)luaL_checkudata(lua, 1, tdigest_metatable);
  double q = luaL_checknumber(lua, 2);
  double x = (*data)->quantile(q);

  lua_pushnumber(lua, x);
  return 1;
}

static int lua_tdigest_tostring(lua_State *lua) {
  int n = lua_gettop(lua);
  luaL_argcheck(lua, n == 1, 0, "incorrect number of arguments");

  TDigest **data = (TDigest **)luaL_checkudata(lua, 1, tdigest_metatable);
  stringstream tmp;
  tmp << **data;
  lua_pushfstring(lua, "%s", tmp.str().c_str());
  return 1;
}

static int lua_tdigest_gc(lua_State *lua) {
  TDigest **data = (TDigest **)luaL_checkudata(lua, 1, tdigest_metatable);
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
  luaL_newmetatable(lua, tdigest_metatable);
  lua_pushvalue(lua, -1);  /* duplicate the metatable */
  lua_setfield(lua, -2, "__index");

  luaL_setfuncs(lua, tdigest_m, 0);
  luaL_newlib(lua, tdigest_f);
  return 1;
}
