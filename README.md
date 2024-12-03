The followings are implemented: 

Geodesic Computation:

- I have implemented 2 methods: Floyd Warshall Algorithm by Arrays and Dijkstra's Algorithm by Heaps
- 
  I measured the execution time for finding all geodesics (from each vertex to all the other vertices).
  
  Heap implementation is much quicker than Floyd Warshall's.
  
  I made it draw the geodesic with a red line. I could not adjust the line thickness since it was giving error (I don't know why), so I need to show it with a thin red line, sorry for the line thickness.. <br/>

  
SDF Segmentation:

Implemeted the paper "[Consistent Mesh Partitioning and Skeletonisation using the Shape Diameter Function.](https://link.springer.com/content/pdf/10.1007/s00371-007-0197-5.pdf)" [1]

- I created some rays from the center of each triangle, but giving a direction to those rays with a specified angle was not an easy task.

  Here is my own direction determining algorithm:

  First I took the corner points and mid-points of edges of the triangle. These give me 6 directions from the center, but on the triangle plane of course.. Then I went from the center in the direction of inner normal as much as 1 unit, and determined the point x. Then I went from the point x in the direction of one of the 6 directions as much as the tangent of the determined angle, and found the target point p. Finally, I found the ray the direction from the center of the triangle to the point p.
  
  First, I sent 1 ray in the direction of inner normal, but this did not give any good results. Next, addition to that ray, I have sent 30 more rays in such a way that taking 6 directions on the triangle plane and 5 angles between 0 and 45 (I chose 9,18,27,36,45). Hence, I generated 5*6=30 different days in different directions in a homegoneous way.
  
  My results are neither too bad, nor good enough. I think that some triangles do not have a consistent sdf value. For instance, in the ceraut model, on the front part of chest there are some triangles with different color than its neighbors. Here, my guess is, the ray in the inner normal direction goes to too far, till the tail of the ceraut whereas the other rays (from the same triangle) goes to much more closer triangles (to the nearby chest triangles). Hence the sdf values are too inconsistent for such a triangle when we compare the other nearby triangles. Moreover, for the hand model, the width from the front to end of the below part is very similar to the width from the front to end of the fingers. Although the widths from the left to right are different for the below part of hand and fingers of hand, this does not reflect to the sdf values. Therefore, I cound not obtain a multi-colored segmentation for the hand model.
  
  By the way, I used histogram clustering to partition the sdf values.. An advantage of SDF segmentation to random walk segmentation is SDF segmentation is much quicker than the other. The disadvantage is results are worse. <br/>

  
RANDOM WALK Segmentation

Implemeted the paper "[Rapid and effective segmentation of 3D models using random walks.](https://cg.cs.tsinghua.edu.cn/papers/cagd_2009_segmentation.pdf)" [2]

- This really gives good results. I have implemented FPS seed selection (it is great), and random seed selection (good).. However, the linear system solution exists for the 3-neighbor triangles. Interestingly, our many models include some triangles having more than 3 neighbors.. Hence the system solution vanishes, and I can not get a result.. Probably I could not get result for some meshes because of the reason that I said above, but still an other possible reason may be  the long computation time of matrix calculations. I have used Eigen's SparseQR sparse matrix solving library. To sum up, I could not get any result for some meshes even though I wait more than 2 hours. I have obtained really good results for the horse.
  
  YOU CAN FIND MY SCREEN SHOTS IN THE "Screenshots" FOLDER.
  
  YOU CAN BUILD MY CODE BY JUST TYPING "make" command.

  YOU CAN RUN MY CODE BY JUST TYPING "./Main" command.
  
  I have implemented in Linux Platform. <br/>

  
ADDITONAL THINGS

- I made it draw rays that going through the cone in the SDF segmentation (The rays are drawn in the reverse direction to be able to be seen).

- Also, you can adjust the maximum angle that the ray makes with the innerNormal inside the cone, and number of angle intervals:

  For example, if maximum angle is 60 and angleCount is 10, then I will send the rays with the following angles: 6,12,18,24,30,36,42,48,54,60

  Another example, if maximum angle is 45 and angleCount is 5, then I will send the rays with the following angles: 9,18,27,36,45

  You can adjust them by changing the values of AngleCount and maxAngle variables in lines 198 and 199, resp., inside Mesh.cpp file.
  (I did not put this feature usable through standard input/output mechanism on Terminal in order not to make the User get bored and confused with many optionalities). 

  In conclusion, I have tried the SDF with different number of rays with different amount of angles, but the results did not change much. <br/>


[1] Shapira, L., Shamir, A., & Cohen-Or, D. (2008). Consistent mesh partitioning and skeletonisation using the shape diameter function. The Visual Computer, 24, 249-259.

[2] Lai, Y. K., Hu, S. M., Martin, R. R., & Rosin, P. L. (2009). Rapid and effective segmentation of 3D models using random walks. Computer Aided Geometric Design, 26(6), 665-679.


  
  
  
  
