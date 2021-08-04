# Jericho Lawson
# 29 July 2021
# Sample Entropy Test

data = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/04015MITexpanded.csv")

sdrr = sd(data$RRLength)
lengths = data$RRLength

# finds standardized SampEnt for segment of RR intervals
sampEntStand = function(x, m = 2, r = 0.1 * sdrr, n = length(x)){
  # solve for A and B
  A = sampEntCheck(x, m + 1, r, n)
  B = sampEntCheck(x, m, r, n)
  
  #return(c(A, B)) # debugging purposes
  
  # return standardized sample entropy
  return(-log(A[3] / B[3]))
}

# finds the number of pairs that are lower than threshold for SampEnt
sampEntCheck = function(x, m, r, n){
  # tracks successful differences
  count = 0
  
  # finds max possible differences
  max = choose(n - m + 1, 2)
  
  # pairs each possible m-pack to another m-pack
  for (i in 1:(n + 1 - m - 1)){
    i_pair = x[i:(i + m - 1)]
    for (j in (i + 1):(n + 1 - m)){
      j_pair = x[j:(j + m - 1)]
      count = count + sum(max(abs(i_pair - j_pair)) <= r)
    }
  }
  
  # returns count, maximum, and rate
  return(c(count, max, count / max))
}

sampEntStand(x = lengths[1:300])
