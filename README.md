# 3D Model Viewer
Using olcPixelGameEngine as a rasterizer, start by reading all the vertex coordinates and geometric connections to create polygons in the loaded .obj file. If it’s textured, read and store all the UV maps in those polygons as well. The program is only capable of rendering the model with one texture file, but UV maps allow for several different textures to be arranged in one image to give the appearance of multiple materials.

First, each polygon is translated by the world offset if there is any, then their normal is calculated. If the dot product of the polygon’s normal and its vector towards the camera is greater than 0, it means its back is facing the camera and should be culled.

Calculate the distance of its midpoint to the camera for lighting purposes, then transform its vertices along the view matrix, which is the camera’s point-at matrix (rotation+translation of a point to another point) but inverted so that the process applies to the polygon instead. The polygon has now transformed from worldspace to viewspace.

If the resulting viewspace polygon clips past the camera, clip off the section that does, otherwise that clipping section will be considered infinite in size and crash the program. Clipping is the process of separating a triangle into two or three smaller triangles depending on whether one or two points lie outside of the clipping plane. The resulting inner triangles are pushed into the rendering pipeline, while the original polygon is discarded.

All the polygons in the rendering pipeline then transform along the projection matrix, making them 2d shapes on the screen. The resulting W value in each vertex now represents their distance from the camera, and thus becomes useful for depth buffering. It’s also critical for aligning the texture’s UV maps correctly, otherwise they would distort with the camera’s perspective like a PS1 game (This problem is called affine texture mapping and is prevalent with older CPU-heavy renderers. and even emulated in some indie games).

The resulting 2d vertex coordinates are normalized into cartesian space, and their X/Y values have to be un-inverted because of the previous projection transformation. These cartesian coordinates are multiplied by the program window’s height and width to finally complete the 2d projection. The engine wipes the display with a background color before rendering the new projection.

For each projected polygon that clips with the side of the screen, clip them against every side of the screen in order of top->bottom->left->right. Discard the original polygon and keep all the new inner triangles. At last, none of the triangles in the rendering pipeline will cause performance problems and can be rendered.

Rendering is done on a pixel-by-pixel basis, starting from the top of the triangle to the bottom, with each row of pixels drawn left to right. This is done by separating the triangle into two halves with a horizontal slice if necessary (a fully horizontal side means one of these “halves” is nonexistent). The vertices are sorted locally in the rendering function such that the aforementioned rendering order is maintained.

The triangle is iterated through with floating point steps calculated beforehand both vertically and horizontally. For each row, the locations of the left and rightmost pixels are determined first, then each pixel on the row is iterated through. The resulting floating point coordinates are converted to integers in screen space (This causes some inaccuracies with specific angles and distances). The corresponding locations of each of these pixels is tested on the polygon’s UV map to grab the correct color from the texture. However, if the current pixel’s calculated depth value is behind that of a previously rendered pixel, skip drawing it. The model is now fully projected on the screen with correctly oriented textures and any pixel out of view culled.


# Maze solver output reading
A maze solver outputs the coordinates of each wall in a generated maze to a .txt file, which the 3D engine reads to build a model of the maze by placing a cube in each location. Once the AI completes the maze, the program outputs the coordinate of every location it has visited in sequential order to another .txt file. After building the model, the 3D engine pushes the sequential AI locations into a vector and reads them in order.

The movement between locations is interpolated over a short time period with a brief pause between each. The camera location is simply set to that of the AI’s, while the orientation at each point in time is generated based on the cardinal direction differential between each step. For example, a movement while facing north from (1,1) towards (2,1) means that the camera should turn right, and that its new orientation is to the east. The brief pause between each movement is utilized to show the camera turning smoothly when one occurs. 180 degree turns are simply ignored and displayed as backwards movement since they are usually a result of backtracking (And it would look like your neck snapped with how jarring the movement would be). This backwards state continues until the AI makes another 180, or takes any left or right turn. The AI will then continually repeat the moving process with the given output of coordinates.


Demo of the maze solver, texturing, and proximity lighting:
https://drive.google.com/file/d/1hBnIPKRqEc05AKIIc8NxgrWzMVloLlOI/view?usp=sharing
