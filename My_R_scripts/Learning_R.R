
######################### cat() and print() ########################

cat("a","b","c")
cat("a","b","c",sep=",",fill=TRUE)
cat("a","b","c",sep="\t",fill=TRUE)

# specify a width after which a new line is inserted
cat("a","b","c",fill=2)
cat("a","b","c",fill=4)




######################### manipulate number #########################

round(3.1415, digits = 2)
# 3.14




#################### for loop and if statements ####################

sum = 0

for (i in 1:5){
  
  if (i >= 3){
    sum = sum + i
  }
  
}

print(sum)




#################### for loop and if statements ####################

# class()      - what kind of object is it (high-level)?
# typeof()     - what is the object’s data type (low-level)?
# length()     - how long is it? What about two dimensional objects?
# attributes() - does it have any metadata?


list_demo = c(1, 2.5, -3, '4', 'good', NA, NaN, TRUE)

for (i in list_demo){
  
  print(i)
  print(class(i))
  cat('\n')
  
}


print(class(1.1))












