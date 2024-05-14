#
# Set overall margins for the combined set of plots and size them
# to generate a requested inter-plot spacing
#
set multiplot layout 2,2
plot sin(x) lt 1
plot cos(x) lt 2
plot sin(2*x) lt 3
plot cos(2*x) lt 4
unset multiplot