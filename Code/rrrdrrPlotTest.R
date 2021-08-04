# Jericho Lawson
# 27 July 2021
# NEC Rate Plot

# libraries
library(ggplot2); #library(plotly)

# function to do plus/minus operation
plusminus = function(x, diff){
  if (length(x) != 2){
    NULL
  } else {
    c(x[1] - diff, x[2] + diff)
  }
}

# MIT data
data2 = read.csv("C:/Users/Mario/Downloads/merged_segmentsNEW.csv")
data3 = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/mit-bih-atrial-fibrillation-database-1.0.0/files/merged_segmentsNEW.csv")
data1 = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/all_trans_dataMIT.csv")
head(data1)

round(data1$Time[4:1000], 2) == round(data3$secs[1:997], 2)
data1[272:277, ]; data3[269:274, ]

# subject's heartbeat data
test = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/04015MITexpanded.csv")
head(test)

# function to determine NEC cells for a subject
# notes: 25 ms, adds equal remaining space to leftmost and rightmost points to get equal 25ms segments
gridplace = function(data = test, varpick = c(3, 9), rows = c(134, 179), int = 0.025){
  xl = plusminus(range(data[rows[1]:rows[2], varpick[1]]), (int - (diff(range(data[rows[1]:rows[2], varpick[1]])) %% int)) / 2)
  yl = plusminus(range(data[rows[1]:rows[2], varpick[2]]), (int - (diff(range(data[rows[1]:rows[2], varpick[2]])) %% int)) / 2)
  
  # grid lines
  xcuts = seq(xl[1], xl[2], by = int); ycuts = seq(yl[1], yl[2], by = int)
  
  # cell counts
  grid = matrix(0, length(ycuts) - 1, length(xcuts) - 1)
  
  # tally up counts of each cell
  for (p in 1:length(data[rows[1]:rows[2], varpick[1]])){
    coords = c(ceiling((data[rows[1]:rows[2], varpick[1]][p] - xl[1]) / int), ceiling((yl[2] - data[rows[1]:rows[2], varpick[2]][p]) / int))
    grid[coords[2], coords[1]] = grid[coords[2], coords[1]] + 1
  }
  
  ## ggplot
  draw = ggplot(data[rows[1]:rows[2], ], aes(x = RRLength, y =  RRDiff)) + geom_point(size = 0.005, col = "red") + 
    theme(panel.grid.minor = element_line(colour = "lightblue", size = 0.005)) + 
    scale_x_continuous(minor_breaks = xcuts) + scale_y_continuous(minor_breaks = ycuts)
  print(draw)
  
  return(grid)
}

# results
a = gridplace()
sum(a != 0)

# this process has to be done when converting heartbeat data into segment data
