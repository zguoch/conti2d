
@mainpage conti2d


# Description

**conti** is a C++ program to calculate upward  and downward continuation of potential field data.
 
![](Geometry_UWC.svg)

# Synopsis

conti2d input.grd -Goutput.grd 

# Required Arguments

* *input.grd*: 待延拓的数据文件
* -**G***output.grd*：延拓结果文件
* -**T**[*h1* | *topo1.grd*]: 输入数据所在平面的高程**h1**，或者所在曲面的高程数据文件**topo1.grd**，必须与输入数据具有相同的网格。
* -**H**[*h2* | *topo2.grd*]: 输出数据所在平面的高程**h2**, 或者所在曲面的高程数据文件**topo2.grd**，同样与输入数据具有相同的网格。

# Optional Arguments¶

* -**E***extNum*: 扩边*extNum*个点距
* -**t***ThreadNum*：指定参与计算的线程数，默认使用所有线程
* -e: 准确解的文件名，必须与输入数据具有相同的网格结构
* -D[+**T***lambda* | +**I***kmax* | +**C***delta* | +**L***kmax*]: 表示向下延拓。**+T*lambda***表示Tikhonov正则化方法，*lambda*指定正则化参数；**I**表示积分迭代方法，*kmax*表示迭代次数；**C**表示使用CGLS方法，*delta*表示收敛条件，即将解上延与输入数据作差后的二范数与输入数据的二范数之比小于某个值; **L**表示Landweber迭代方法，*kmax*表示迭代次数。
* -**f**：表示基于FFT在频率域计算，默认在空间域计算。但是频率域计算只限平面对平面的延拓。
* 