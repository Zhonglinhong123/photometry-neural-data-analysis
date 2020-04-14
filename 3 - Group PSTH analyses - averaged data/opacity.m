% calculates translucent rgb values for given rgb input + o% opacity

function [opacity] = opacity (o, r, g, b)

r_o = 255;
g_o = 255;
b_o = 255;

opacity = [((1-o)*r + o*r_o)/255, ((1-o)*g + o*g_o)/255, ((1-o)*b + o*b_o)/255];