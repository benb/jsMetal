#!/usr/bin/env ruby
require 'chunky_png'
png = ChunkyPNG::Image.new(100, 10, ChunkyPNG::Color('gray'))
png.save('0.png')
(0...100).each do |i| 
        (0...10).each do |j|
                png[i,j]=ChunkyPNG::Color('red')
        end
        png.save((i+1).to_s + '.png')
end

