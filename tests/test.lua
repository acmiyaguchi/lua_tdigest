local tdigest = require "tdigest"
local math = require "math"

local function box_muller(scale, base)
    return base + scale*math.sqrt(-2 * math.log(math.random())) * math.cos(2 * math.pi * math.random()) / 2
end

local function quantile(sample, q)
   return sample[math.floor(q*#sample)]
end

local td = tdigest.new()
local sample = {}

for i = 1, 10000 do
   local x = box_muller(5, 2) + box_muller(15, 4)
   sample[#sample + 1] = x
   td:add(x)
end

table.sort(sample);
for i = 0.1, 0.9, 0.1 do
   assert(math.abs(quantile(sample, 0.5) - td:quantile(0.5)) < 0.01)
end
