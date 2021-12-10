% Wrap.m
%
% Aaron James Lemmer
% November 7, 2013

function [wrapped_phase_difference] = Wrap(phase_difference)

wrapped_phase_difference = atan2(sin(phase_difference), cos(phase_difference));