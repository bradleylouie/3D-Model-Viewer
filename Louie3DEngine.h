#ifndef HEADER_8D040547EFE6B8A3
#define HEADER_8D040547EFE6B8A3

#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include <fstream>
#include <sstream>
#include <algorithm> //for painter's algorithm which should be replaced with depth buffers later
#include <math.h>
using namespace std;

struct vec2d
{
	float u = 0;
	float v = 0;
	float w = 1;
	vec2d () {}

	vec2d ( float uval, float vval )
	{
		u = uval; v = vval;
	}
};

struct vec3d
{
	float x = 0;
	float y = 0;
	float z = 0;
	float w = 1; //used for 4x4 matrix calculations
    //doubles are unnecessarily precise, so floats will suffice

    vec3d () {}

    vec3d (float xval, float yval, float zval)
    {
        x = xval; y = yval; z = zval;
    }

	vec3d operator + (const vec3d& rhs)
	{
		return { x + rhs.x, y + rhs.y, z + rhs.z };
	}
		void operator += (const vec3d& rhs)
		{
			x += rhs.x;
			y += rhs.y;
			z += rhs.z;
		}

	vec3d operator - (const vec3d& rhs)
	{
		return { x - rhs.x, y - rhs.y, z - rhs.z };
	}
		void operator -= (const vec3d& rhs)
		{
			x -= rhs.x;
			y -= rhs.y;
			z -= rhs.z;
		}

	vec3d operator * (const float& scalar)
	{
		return { x * scalar, y * scalar, z * scalar };
	}

	vec3d operator / (const float& scalar)
	{
		return { x / scalar, y / scalar, z / scalar };
	}

	bool operator == (const vec3d& rhs)
	{
	    if (x == rhs.x && y == rhs.y && z == rhs.z && w == rhs.w)
            return true;
        return false;
	}
        bool operator != (const vec3d& rhs)
        {
            if (x == rhs.x && y == rhs.y && z == rhs.z && w == rhs.w)
                return false;
            return true;
        }

	float DotProduct(const vec3d& rhs)
	{
		return x*rhs.x + y*rhs.y + z*rhs.z;
	}

	vec3d CrossProduct(const vec3d& rhs)
	{
		return { y*rhs.z - z*rhs.y, z*rhs.x - x*rhs.z, x*rhs.y - y*rhs.x };
	}

	float GetLength()
	{
		return sqrtf(this->DotProduct(*this));
	}

	vec3d Normalize()
	{
		float l = this->GetLength();
		return { x / l, y / l, z / l };
	}

};

//given two points, a plane's position, and its normal, return the position where they intersect
//also set the value of t to the normalized value of the intersection point along the vector, e.g. in the middle = 0.5
vec3d VectorIntersectPlane(vec3d& plane_pos, vec3d& plane_n, vec3d& line_start, vec3d& line_end, float& t)
{
	//Ensure the normal is normalized
	vec3d plane_normal = plane_n.Normalize();
	float plane_d = -plane_normal.DotProduct(plane_pos);
	float ad = line_start.DotProduct(plane_normal);
	float bd = line_end.DotProduct(plane_normal);
	t = (-plane_d - ad) / (bd - ad);
	vec3d line_length = line_end - line_start;
	return line_start + (line_length * t);
}

struct poly
{
	vec3d vertex[3]; //hopefully in clockwise order
	vec2d texture_vertex[3];
	float luminance = 1.0f;
	olc::Pixel color;

	//Given a plane's position and its normal, return up to two new polys based on how this poly clips with the plane
	int ClipAgainstPlane(vec3d plane_pos, vec3d plane_n, poly& out_tri1, poly& out_tri2)
	{
		//Ensure the normal is normalized
		vec3d plane_normal = plane_n.Normalize();

		//Returns shortest signed distance from given point to the plane
		auto dist = [&](vec3d &vertex)
		{
			return (plane_normal.x*vertex.x + plane_normal.y*vertex.y + plane_normal.z*vertex.z - plane_normal.DotProduct(plane_pos));
		};

		//Group points by whether they lie "beyond" or "before" the plane
		vec3d* inside_points[3];  int inside_point_count = 0;
		vec3d* outside_points[3]; int outside_point_count = 0;
		vec2d* inside_tex[3];
		vec2d* outside_tex[3];

		//Store signed distance of each point in the poly to the plane
		float d0 = dist(vertex[0]);
		float d1 = dist(vertex[1]);
		float d2 = dist(vertex[2]);

		//Count and store each point based on whether they lie "beyond" the plane
		if (d0 >= 0) {
            inside_points[inside_point_count] = &vertex[0];
			inside_tex[inside_point_count] = &texture_vertex[0];
            inside_point_count++;
        } else {
            outside_points[outside_point_count] = &vertex[0];
			outside_tex[outside_point_count] = &texture_vertex[0];
            outside_point_count++;
        }
		if (d1 >= 0) {
            inside_points[inside_point_count] = &vertex[1];
			inside_tex[inside_point_count] = &texture_vertex[1];
            inside_point_count++;
        } else {
            outside_points[outside_point_count] = &vertex[1];
			outside_tex[outside_point_count] = &texture_vertex[1];
            outside_point_count++;
        }
		if (d2 >= 0) {
            inside_points[inside_point_count] = &vertex[2];
			inside_tex[inside_point_count] = &texture_vertex[2];
            inside_point_count++;
        } else {
            outside_points[outside_point_count] = &vertex[2];
			outside_tex[outside_point_count] = &texture_vertex[2];
            outside_point_count++;
        }

		//Cull or split poly based on four possible cases (0, 1, 2, or all points lie within the viewing frustrum)
		if (inside_point_count == 0) //ALL points outside plane, cull the whole poly
		{
			return 0; //Do not use any returned polys
		}

		if (inside_point_count == 3) //ALL points inside plane, do nothing
		{
			out_tri1 = *this;
			return 1; //Use the original poly untouched
		}

		if (inside_point_count == 1 && outside_point_count == 2) //Two points outside plane, clip to smaller triangle
		{
			//Copy appearance info
			//out_tri1.color = olc::RED;
			out_tri1.color = color;
			out_tri1.luminance = luminance;

			//Keep the one inner point
			out_tri1.vertex[0] = *inside_points[0];
			out_tri1.texture_vertex[0] = *inside_tex[0];

			//Form a triangle with the inside point and the intersections of its vectors to the two outside points
			float t; //normal value of the intersection point along the vector
			out_tri1.vertex[1] = VectorIntersectPlane(plane_pos, plane_normal, *inside_points[0], *outside_points[0], t);
			out_tri1.texture_vertex[1].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.texture_vertex[1].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.texture_vertex[1].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;

			out_tri1.vertex[2] = VectorIntersectPlane(plane_pos, plane_normal, *inside_points[0], *outside_points[1], t);
			out_tri1.texture_vertex[2].u = t * (outside_tex[1]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.texture_vertex[2].v = t * (outside_tex[1]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.texture_vertex[2].w = t * (outside_tex[1]->w - inside_tex[0]->w) + inside_tex[0]->w;

			return 1; //Use the new poly
		}

		if (inside_point_count == 2 && outside_point_count == 1) //One point outside plane, clip into TWO smaller triangles
		{
			//out_tri1.color = olc::BLUE;
			//out_tri2.color = olc::GREEN;
            out_tri1.color = color;
            out_tri2.color = color;
            out_tri1.luminance = luminance;
            out_tri2.luminance = luminance;
			//Form two triangles to create the new quad
			float t;
			out_tri1.vertex[0] = *inside_points[0];
			out_tri1.vertex[1] = *inside_points[1];
			out_tri1.texture_vertex[0] = *inside_tex[0];
			out_tri1.texture_vertex[1] = *inside_tex[1];
			out_tri1.vertex[2] = VectorIntersectPlane(plane_pos, plane_normal, *inside_points[0], *outside_points[0], t);
			out_tri1.texture_vertex[2].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.texture_vertex[2].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.texture_vertex[2].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;

			out_tri2.vertex[0] = *inside_points[1];
			out_tri2.vertex[1] = out_tri1.vertex[2];
			out_tri2.texture_vertex[0] = *inside_tex[1];
			out_tri2.texture_vertex[1] = out_tri1.texture_vertex[2];
			out_tri2.vertex[2] = VectorIntersectPlane(plane_pos, plane_normal, *inside_points[1], *outside_points[0], t);
			out_tri2.texture_vertex[2].u = t * (outside_tex[0]->u - inside_tex[1]->u) + inside_tex[1]->u;
			out_tri2.texture_vertex[2].v = t * (outside_tex[0]->v - inside_tex[1]->v) + inside_tex[1]->v;
			out_tri2.texture_vertex[2].w = t * (outside_tex[0]->w - inside_tex[1]->w) + inside_tex[1]->w;

			return 2; //Use the two new polys
		}
		return 0; //never should reach this point
	}
};

struct model
{
    vector<poly> polys;
	bool textured = false;
	bool AIPath = false;

    bool LoadFromObj(string filepath)
    {
        ifstream obj_file;
        obj_file.open(filepath);

        if (!obj_file.is_open())
            return false;

        //stores vertices listed by file
        vector<vec3d> vertices;
        string line;
        vec2d t1 = vec2d(0.0f, 0.0f);
        vec2d t2 = vec2d(1.0f, 0.0f);
        vec2d t3 = vec2d(0.0f, 1.0f);

        while (getline(obj_file, line))
        {
            stringstream ss(line);

            //lines usually start with '#' (comment), 'v', 'f', or 's' (some kind of flag); anything not 'v' or 'f' is considered junk and ignored
            char vertex_or_face;
            ss >> vertex_or_face;

            if (vertex_or_face == 'v') //line represents vertex
            {
				vec3d v;
				ss >> v.x >> v.y >> v.z;
				//e.g. line is 'v 1.0 -1.0 0.0' -> create a vertex at coordinates (1.0, -1.0, 0.0)
				vertices.push_back(v);
            }

			if (vertex_or_face == 'f')  //line represents face
			{
				int c[3]; //contains the connectivity of vertices in clockwise order

				ss >> c[0] >> c[1] >> c[2];
				//e.g. line is 'f 1 2 3' -> create a face from vertices 1, 2, and 3 (but indexed starting from 0)
				//obj's reuse vertices several times to be efficient, as opposed to a soup of polygons
				polys.push_back( {vertices[c[0] - 1], vertices[c[1] - 1], vertices[c[2] - 1], t1, t2, t3} );
			}
        }
        return true;
    }

	bool LoadFromObjTextured(string filepath)
    {
        ifstream obj_file;
        obj_file.open(filepath);

        if (!obj_file.is_open())
            return false;

		textured = true;

        //stores vertices listed by file
        vector<vec3d> vertices;
		vector<vec2d> tex_verts;
        string line;

        while (getline(obj_file, line))
        {
            stringstream ss(line);

            //lines usually start with '#' (comment), 'v', 'f', or 's' (some kind of flag); anything not 'v' or 'f' is considered junk and ignored
            string vertex_or_face;
            ss >> vertex_or_face;

            if (vertex_or_face == "v") //line represents vertex
            {
				vec3d v;
				ss >> v.x >> v.y >> v.z;
				//e.g. line is 'v 1.0 -1.0 0.0' -> create a vertex at coordinates (1.0, -1.0, 0.0)
				vertices.push_back(v);
            }

			if (vertex_or_face == "vt") //line represents texture vertex
            {
				vec2d v;
				ss >> v.u >> v.v;
				//e.g. line is 'vt 0.0 1.0' -> create a texture vertex at (0.0, 1.0)
				tex_verts.push_back(v);
            }

			if (vertex_or_face == "f")  //line represents face
			{
				int c[3]; //contains the connectivity of 3d vertices in clockwise order
				int tc[3]; //contains the connectivity of texture vertices in clockwise order

				//process out the '/' in between the connective values
				string decode[3];
				ss >> decode[0] >> decode[1] >> decode[2];
				for (int i; i < 3; i++)
					{c[i] = decode[i][0]; tc[i] = decode[i][2];}

				//e.g. line is 'f 1/3 2/2 3/1' -> create a face from vertices 1, 2, and 3, using
				//                                texture vertices 3, 2, and 1 (but indexed starting from 0)
				polys.push_back( {vertices[c[0] - 1], vertices[c[1] - 1], vertices[c[2] - 1],
							tex_verts[tc[0] - 1], tex_verts[tc[1] - 1], tex_verts[tc[2] - 1]} );
			}
        }
        return true;
    }

	bool ConstructAIMaze(string filepath)
    {
		static vector<vec2d> t_verts = { vec2d( 0.0f, 0.0f ), vec2d( 1.0f, 0.0f ),
								 		 vec2d( 0.0f, 1.0f ), vec2d( 1.0f, 1.0f ) };
        ifstream maze_file;
        maze_file.open(filepath);

        if (!maze_file.is_open())
            return false;

		AIPath = true;
        //stores vertices listed by file
        string line;

        while (getline(maze_file, line))
        {
            stringstream ss(line);
			float x_pos, z_pos;

			//set coordinate of bottom-center of cube
			ss >> x_pos >> z_pos;
			vector<vec3d> vertices{ vec3d( x_pos-0.5f, 1.0f, z_pos+0.5f ), vec3d( x_pos+0.5f, 1.0f, z_pos+0.5f ),
                                    vec3d( x_pos-0.5f, 1.0f, z_pos-0.5f ), vec3d( x_pos+0.5f, 1.0f, z_pos-0.5f ),
                                    vec3d( x_pos-0.5f, 0.0f, z_pos+0.5f ), vec3d( x_pos+0.5f, 0.0f, z_pos+0.5f ),
                                    vec3d( x_pos-0.5f, 0.0f, z_pos-0.5f ), vec3d( x_pos+0.5f, 0.0f, z_pos-0.5f ) };

			//front face
			polys.push_back({vertices[6], vertices[2], vertices[3], t_verts[2], t_verts[0], t_verts[1]});
			polys.push_back({vertices[6], vertices[3], vertices[7], t_verts[2], t_verts[1], t_verts[3]});
			//back face
			polys.push_back({vertices[5], vertices[1], vertices[0], t_verts[2], t_verts[0], t_verts[1]});
			polys.push_back({vertices[5], vertices[0], vertices[4], t_verts[2], t_verts[1], t_verts[3]});
			//left face
			polys.push_back({vertices[4], vertices[0], vertices[2], t_verts[2], t_verts[0], t_verts[1]});
			polys.push_back({vertices[4], vertices[2], vertices[6], t_verts[2], t_verts[1], t_verts[3]});
			//right face
			polys.push_back({vertices[7], vertices[3], vertices[1], t_verts[2], t_verts[0], t_verts[1]});
			polys.push_back({vertices[7], vertices[1], vertices[5], t_verts[2], t_verts[1], t_verts[3]});
			//bottom face
			polys.push_back({vertices[7], vertices[5], vertices[4], t_verts[2], t_verts[0], t_verts[1]});
			polys.push_back({vertices[7], vertices[4], vertices[6], t_verts[2], t_verts[1], t_verts[3]});
			//top face
			polys.push_back({vertices[2], vertices[0], vertices[1], t_verts[2], t_verts[0], t_verts[1]});
			polys.push_back({vertices[2], vertices[1], vertices[3], t_verts[2], t_verts[1], t_verts[3]});

        }
        return true;
    }
};


vec3d NORTH = vec3d(0,0,1), WEST = vec3d(1,0,0), SOUTH = vec3d(0,0,-1), EAST = vec3d(-1,0,0);
const int NONE = 0, LEFT = -1, RIGHT = 1, BACKTRACK = 2;
const float MOVE_PERIOD = 0.5f;
const float MOVE_INTERVAL = 0.4f;

struct AIPath
{

    int index = 0;
    float step_time = 0;

    vector<vec3d> location;  //location of AI at time t
    vector<float> direction; //cardinal direction in radians of location t towards location t+1; vector size is 1 less than location's
    vector<int> turns;       //whether AI turns 90 degrees left, right, or not at all at time t
    bool ConstructAIPath(string filepath)
    {
        ifstream path_file;
        path_file.open(filepath);

        if (!path_file.is_open())
            return false;

        bool backtracking = false;
        turns.push_back(NONE);
        string line;
        while (getline(path_file, line))
        {
            stringstream ss(line);
			float x_pos, z_pos;

			ss >> x_pos >> z_pos;
            location.push_back({x_pos, 0.5f, z_pos});

            //get 3d vector from previous pos to current pos, then set turn if direction is different from previous
            if (location.size() > 1)
            {
                vec3d dir = location.back() - location.at(location.size()-2);
                if (!backtracking) {
                    if      (dir == NORTH) { direction.push_back( 0.0f );     }
                    else if (dir == EAST ) { direction.push_back( 1.5708f );  }
                    else if (dir == SOUTH) { direction.push_back( 3.14159f ); }
                    else if (dir == WEST ) { direction.push_back( -1.5708f ); }
                    else { runtime_error("Noncardinal direction at " + to_string(location.size()));}}
                else {
                    if      (dir == NORTH) { direction.push_back( 3.14159f ); }
                    else if (dir == EAST ) { direction.push_back( -1.5708f ); }
                    else if (dir == SOUTH) { direction.push_back( 0.0f );     }
                    else if (dir == WEST ) { direction.push_back( 1.5708f );  }
                    else { runtime_error("Noncardinal direction at " + to_string(location.size()));}}

                if (location.size() > 2)
                {
                    vec3d point1 = location.at(location.size()-2) - location.at(location.size()-3);
                    vec3d point2 = location.back() - location.at(location.size()-2);
                    if (point1 != point2)
                    {
                        int turn_dir = GetTurnDirection(point1, point2);
                        turns.push_back(turn_dir);
                        if (turn_dir == BACKTRACK && !backtracking) {
                            backtracking = true;
                            direction.back() += 3.14159;
                        } else if (backtracking) {
                            backtracking = false;
                            turns.back() *= -1;
                            direction.back() += 3.14159;
                        }
                    } else
                        turns.push_back(NONE);
                }
            }
        }
        turns.push_back(NONE);
        return true;
    }

    int GetTurnDirection(vec3d point1, vec3d point2)
    {
        if (point1==point2)      { return NONE; }
        if (point1 == NORTH){
            if (point2 == WEST)  { return LEFT;  }
            if (point2 == EAST)  { return RIGHT; }}
        if (point1 == EAST){
            if (point2 == NORTH) { return LEFT;  }
            if (point2 == SOUTH) { return RIGHT; }}
        if (point1 == SOUTH){
            if (point2 == EAST)  { return LEFT;  }
            if (point2 == WEST)  { return RIGHT; }}
        if (point1 == WEST){
            if (point2 == SOUTH) { return LEFT;  }
            if (point2 == NORTH) { return RIGHT; }}
        return BACKTRACK;
    }

    void UpdateVars(float total_elapsed_time)
    {
        index = (int) ((total_elapsed_time + MOVE_INTERVAL) / (MOVE_PERIOD + MOVE_INTERVAL));
        step_time = fmod(total_elapsed_time, (MOVE_PERIOD + MOVE_INTERVAL));
        if (index >= location.size()){
            index = location.size()-1;
            step_time = 0;
        }
    }

    float GetAngle()
    {
        float t = step_time - MOVE_PERIOD;
        if (index < 2)
            return direction.at(0);
        if ( t > 0 )
        {
            if (turns.at(index-1) == LEFT)
                return direction.at(index-2) - 1.5708 * (t / MOVE_INTERVAL);
            else if (turns.at(index-1) == RIGHT)
                return direction.at(index-2) + 1.5708 * (t / MOVE_INTERVAL);
            else
                return direction.at(index-1);
        }
        return direction.at(index-1);
    }

    vec3d GetPosition()
    {
        if (index == 0)
            return location.at(0);
        if (index == location.size()-1)
            return location.at(index);

        if ( step_time < MOVE_PERIOD )
            return location.at(index-1) + ((location.at(index) - location.at(index-1)) * (step_time/MOVE_PERIOD));
        else
            return location.at(index-1);
    }
};

const int IDENTITY = 0, ROTATEX = 1, ROTATEY = 2, ROTATEZ = 3, TRANSLATION = 4, PROJECTION = 5, POINTAT = 6;
struct mat44f
{
	public:
        float m[4][4] = {0};

        struct Proxy {
            float* _array;

            Proxy(float* _array) : _array(_array) {}

            float operator[](int index) {
                return _array[index];
            }
        };

        Proxy operator []( int i ) { return Proxy(m[i]); }

		mat44f(){}

		//Generate Identity matrix
		mat44f(int type_of_matrix)
		{
			if (type_of_matrix == IDENTITY){
				m[0][0] = 1.0f;
				m[1][1] = 1.0f;
				m[2][2] = 1.0f;
				m[3][3] = 1.0f;
			} else {
				runtime_error("Wrong matrix type or parameters given!");
			}
		}

		//Generate RotateX/Y/Z matrix
		mat44f(int type_of_matrix, const float& angle_rad)
		{
			if (type_of_matrix == ROTATEX) {
				m[0][0] = 1.0f;
				m[1][1] = cosf(angle_rad);
				m[1][2] = sinf(angle_rad);
				m[2][1] = -sinf(angle_rad);
				m[2][2] = cosf(angle_rad);
				m[3][3] = 1.0f;
			} else if (type_of_matrix == ROTATEY) {
				m[0][0] = cosf(angle_rad);
				m[0][2] = sinf(angle_rad);
				m[2][0] = -sinf(angle_rad);
				m[1][1] = 1.0f;
				m[2][2] = cosf(angle_rad);
				m[3][3] = 1.0f;
			} else if (type_of_matrix == ROTATEZ) {
				m[0][0] = cosf(angle_rad);
				m[0][1] = sinf(angle_rad);
				m[1][0] = -sinf(angle_rad);
				m[1][1] = cosf(angle_rad);
				m[2][2] = 1.0f;
				m[3][3] = 1.0f;
			} else {
				runtime_error("Wrong matrix type or parameters given!");
			}
		}

		//Generate Translation matrix
		mat44f(int type_of_matrix, const float& x, const float& y, const float& z)
		{
			if (type_of_matrix == TRANSLATION) {
				*this = mat44f(IDENTITY);
				m[3][0] = x;
				m[3][1] = y;
				m[3][2] = z;
			} else {
				runtime_error("Wrong matrix type or parameters given!");
			}
		}

		//Generate Projection matrix
		mat44f(int type_of_matrix, const float& fov_deg, const float& aspect_ratio, const float& znear, const float& zfar)
		{
			if (type_of_matrix == PROJECTION) {
				float inverse_fov_func = 1.0f / tanf(fov_deg * 0.5f * 0.0174533f);
				m[0][0] = aspect_ratio * inverse_fov_func;
				m[1][1] = inverse_fov_func;
				m[2][2] = zfar / (zfar - znear);
				m[3][2] = (-zfar * znear) / (zfar - znear);
				m[2][3] = 1.0f;
				m[3][3] = 0.0f;
			} else {
				runtime_error("Wrong matrix type or parameters given!");
			}
		}

		//Generate Pointat matrix
		//For creating a rotation+translation matrix mainly for the camera
		//Can also be used to rotate objects with forward vectors towards a location
		mat44f(int type_of_matrix, vec3d& pos, vec3d& target, vec3d& up)
		{
			if (type_of_matrix == POINTAT) {
				vec3d new_forward = target - pos;
				new_forward = new_forward.Normalize();

				vec3d a = new_forward * up.DotProduct(new_forward);
				vec3d new_up = up - a;
				new_up = new_up.Normalize();

				vec3d new_right = new_up.CrossProduct(new_forward);

				m[0][0] = new_right.x;	    m[0][1] = new_right.y;  	m[0][2] = new_right.z;	    m[0][3] = 0.0f;
				m[1][0] = new_up.x;		    m[1][1] = new_up.y;		    m[1][2] = new_up.z;		    m[1][3] = 0.0f;
				m[2][0] = new_forward.x;	m[2][1] = new_forward.y;	m[2][2] = new_forward.z;	m[2][3] = 0.0f;
				m[3][0] = pos.x;			m[3][1] = pos.y;			m[3][2] = pos.z;			m[3][3] = 1.0f;
			} else {
				runtime_error("Wrong matrix type or parameters given!");
			}
		}

		//matrix multiply matrix
		mat44f operator * (mat44f& rhs)
		{
			mat44f result = mat44f();
			for (int r = 0; r < 4; r++)
				for (int c = 0; c < 4; c++)
					result.m[r][c] = m[r][0]*rhs[0][c] + m[r][1]*rhs[1][c] + m[r][2]*rhs[2][c] + m[r][3]*rhs[3][c];
			return result;
		}

		//matrix multiply 3d vector
		vec3d operator * (const vec3d& v)
		{
			vec3d result;
			result.x = v.x*m[0][0] + v.y*m[1][0] + v.z*m[2][0] + v.w*m[3][0];
			result.y = v.x*m[0][1] + v.y*m[1][1] + v.z*m[2][1] + v.w*m[3][1];
			result.z = v.x*m[0][2] + v.y*m[1][2] + v.z*m[2][2] + v.w*m[3][2];
			result.w = v.x*m[0][3] + v.y*m[1][3] + v.z*m[2][3] + v.w*m[3][3];
			return result;
		}

		//@PRE: rotation and/or translation matrices only
		//@POST: result = inverse of this
		//Mainly used with camera's pointat matrix to invert the functionality so that all objects will transform towards the camera
		//This inverted pointat matrix becomes known as a lookat matrix; objects in worldspace multiply by this matrix to transform into viewspace
		mat44f Invert()
		{
			mat44f result = mat44f();
			result.m[0][0] = m[0][0]; 	result.m[0][1] = m[1][0]; 	result.m[0][2] = m[2][0]; 	result.m[0][3] = 0.0f;
			result.m[1][0] = m[0][1]; 	result.m[1][1] = m[1][1]; 	result.m[1][2] = m[2][1]; 	result.m[1][3] = 0.0f;
			result.m[2][0] = m[0][2]; 	result.m[2][1] = m[1][2]; 	result.m[2][2] = m[2][2]; 	result.m[2][3] = 0.0f;

			result.m[3][0] = -(m[3][0] * result[0][0] + m[3][1] * result[1][0] + m[3][2] * result[2][0]);
			result.m[3][1] = -(m[3][0] * result[0][1] + m[3][1] * result[1][1] + m[3][2] * result[2][1]);
			result.m[3][2] = -(m[3][0] * result[0][2] + m[3][1] * result[1][2] + m[3][2] * result[2][2]);
			result.m[3][3] = 1.0f;
			return result; //treat this as a blackbox i guess...
		}
};

class Louie3DEngine : public olc::PixelGameEngine
{
public:
	Louie3DEngine()
	{
		sAppName = "3D Demo";
	}


private:
	model loaded_mesh;
	AIPath solver_path;
	olc::Sprite texture_spr;
	mat44f projection_matrix;	//converts viewspace to screenspace
	vec3d cam_location;	        //location of camera in the world
	vec3d look_direction;	    //vector of the camera's looking direction
	vec3d right_direction;
	float yaw;		            //sideways camera rotation
	float pitch;	            //Not used yet, for implementing mouselook
	float fov = 90.0f;
	float znear = 0.1f;
	float zfar = 1000.0f;
	float t = 0.0f;
	olc::Pixel bg_color = olc::BLACK;
	olc::Pixel fl_color = olc::BLACK;
	vec3d light_source;
	float current_height = 0.9f; //starting height for AI maze traversing
	float* depth_buffer;         //pseudo-2D array containing depth values for every pixel on the display
	float view_dist_inv;

public:
	bool OnUserCreate() override
	{
		loaded_mesh.ConstructAIMaze("MazeOut.txt");
		solver_path.ConstructAIPath("AgentHistory.txt");
		texture_spr.LoadFromFile("brick.png");
		view_dist_inv = 0.25f;

        //view_dist_inv = 0.012f;
        //texture_spr.LoadFromFile("snow.png");
		//loaded_mesh.LoadFromObj("mountains.obj");

        loaded_mesh.textured = true;

		projection_matrix = mat44f(PROJECTION, fov, (float)ScreenHeight() / (float)ScreenWidth(), znear, zfar);
		//light_source = light_source.Normalize();
		depth_buffer = new float[ScreenWidth()*ScreenHeight()];

		return true;
	}

	bool OnUserUpdate(float elapsed_time) override
	{
        float velocity = 4.0f * elapsed_time;
        mat44f camera_yaw;
		if (loaded_mesh.AIPath)
		{
			solver_path.UpdateVars(t);
			t += elapsed_time;
			cam_location = solver_path.GetPosition();
			cam_location.y += current_height;

			if (GetKey(olc::Key::SPACE).bHeld)
				current_height += velocity;	//Ascend

			if (GetKey(olc::Key::CTRL).bHeld)
				current_height -= velocity;	//Descend

			camera_yaw = mat44f(ROTATEY, solver_path.GetAngle());
		}
		else
		{
			float vertical_velocity = 4.0f * elapsed_time;
			vec3d forward_velocity = look_direction * velocity;
			vec3d right_velocity = right_direction * velocity;
			float turn_rate = 2.0f * elapsed_time;

			if (GetKey(olc::Key::SPACE).bHeld)
				cam_location.y += vertical_velocity; //Ascend

			if (GetKey(olc::Key::CTRL).bHeld)
				cam_location.y -= vertical_velocity; //Descend

			if (GetKey(olc::Key::W).bHeld)
				cam_location += forward_velocity; //Move forward

			if (GetKey(olc::Key::S).bHeld)
				cam_location -= forward_velocity; //Move backward

			if (GetKey(olc::Key::LEFT).bHeld)
				cam_location -= right_velocity;	//Strafe left

			if (GetKey(olc::Key::RIGHT).bHeld)
				cam_location += right_velocity;	//Strafe right

			if (GetKey(olc::Key::A).bHeld)
				yaw -= turn_rate; 				//Turn left

			if (GetKey(olc::Key::D).bHeld)
				yaw += turn_rate;				//Turn right

			camera_yaw = mat44f(ROTATEY, yaw);
		}

		mat44f world_translation = mat44f(TRANSLATION, 0.0f, 0.0f, 0.0f);
		mat44f world_matrix = mat44f(IDENTITY) * world_translation ; //form the translated world matrix
		//Camera setup
		vec3d cam_up = { 0,1,0 };
		vec3d cam_target = { 0,0,1 };

		mat44f camera_yaw_right = mat44f(ROTATEY, yaw+1.5708f);
		look_direction = camera_yaw * cam_target;
		right_direction = camera_yaw_right * cam_target;
		cam_target = cam_location + look_direction;
		light_source = cam_location;
		mat44f camera_matrix = mat44f(POINTAT, cam_location, cam_target, cam_up);

		mat44f view_matrix = camera_matrix.Invert();

		vector<poly> polys_to_rasterize;
		for (auto tri : loaded_mesh.polys)
		{
			poly tri_projected, tri_transformed, tri_viewed;

			//transform every triangle to world space
			tri_transformed.vertex[0] = world_matrix * tri.vertex[0];
			tri_transformed.vertex[1] = world_matrix * tri.vertex[1];
			tri_transformed.vertex[2] = world_matrix * tri.vertex[2];
			tri_transformed.texture_vertex[0] = tri.texture_vertex[0];
			tri_transformed.texture_vertex[1] = tri.texture_vertex[1];
			tri_transformed.texture_vertex[2] = tri.texture_vertex[2];

			//calculate the normals of each triangle
			vec3d line1 = tri_transformed.vertex[1] - tri_transformed.vertex[0];
			vec3d line2 = tri_transformed.vertex[2] - tri_transformed.vertex[0];
			vec3d normal = line1.CrossProduct(line2);
			normal = normal.Normalize();

			vec3d vec_poly_to_camera = tri_transformed.vertex[0] - cam_location;

			//only render if poly's normal is facing the camera
			if (normal.DotProduct(vec_poly_to_camera) < 0.0f)
			{
				//get alignment of surface normal to the light source (0.1 ~ 1.0)
				vec3d midpoint = (tri_transformed.vertex[0]+tri_transformed.vertex[1]+tri_transformed.vertex[2])/3.0f;
				float dist_from_source = sqrt((light_source.x - midpoint.x)*(light_source.x - midpoint.x)
                                              +(light_source.y - midpoint.y)*(light_source.y - midpoint.y)
                                              +(light_source.z - midpoint.z)*(light_source.z - midpoint.z));
				float lum = min(1.0f, max(1.0f - dist_from_source*view_dist_inv, 0.033f));
				tri_transformed.luminance = lum;
                tri_transformed.color = olc::Pixel(255.0f*lum, 255.0f*lum, 255.0f*lum);

				//transform the world space triangles to view space
				tri_viewed.vertex[0] = view_matrix * tri_transformed.vertex[0];
				tri_viewed.vertex[1] = view_matrix * tri_transformed.vertex[1];
				tri_viewed.vertex[2] = view_matrix * tri_transformed.vertex[2];
				tri_viewed.texture_vertex[0] = tri_transformed.texture_vertex[0];
				tri_viewed.texture_vertex[1] = tri_transformed.texture_vertex[1];
				tri_viewed.texture_vertex[2] = tri_transformed.texture_vertex[2];
				tri_viewed.color = tri_transformed.color;
				tri_viewed.luminance = tri_transformed.luminance;

				//clip against znear plane to prevent infinite sizes
				poly clipped[2];
				int num_clipped_polys = tri_viewed.ClipAgainstPlane({ 0.0f, 0.0f, znear }, { 0.0f, 0.0f, 1.0f }, clipped[0], clipped[1]);

				for (int n = 0; n < num_clipped_polys; n++) //n==1 if there was no clipping
				{
					//project from 3d to 2d space; clipped[0] is the original poly if there was no clipping
					tri_projected.vertex[0] = projection_matrix * clipped[n].vertex[0];
					tri_projected.vertex[1] = projection_matrix * clipped[n].vertex[1];
					tri_projected.vertex[2] = projection_matrix * clipped[n].vertex[2];
					tri_projected.texture_vertex[0] = clipped[n].texture_vertex[0];
					tri_projected.texture_vertex[1] = clipped[n].texture_vertex[1];
					tri_projected.texture_vertex[2] = clipped[n].texture_vertex[2];
					tri_projected.color = clipped[n].color;
					tri_projected.luminance = clipped[n].luminance;


                    //use the W value of the vertex to scale the texture coordinate into perspective before normalizing the vertex
                    tri_projected.texture_vertex[0].u /= tri_projected.vertex[0].w;
					tri_projected.texture_vertex[1].u /= tri_projected.vertex[1].w;
					tri_projected.texture_vertex[2].u /= tri_projected.vertex[2].w;

					tri_projected.texture_vertex[0].v /= tri_projected.vertex[0].w;
					tri_projected.texture_vertex[1].v /= tri_projected.vertex[1].w;
					tri_projected.texture_vertex[2].v /= tri_projected.vertex[2].w;

					tri_projected.texture_vertex[0].w /= tri_projected.vertex[0].w;
					tri_projected.texture_vertex[1].w /= tri_projected.vertex[1].w;
					tri_projected.texture_vertex[2].w /= tri_projected.vertex[2].w;

					//normalize into cartesian space
					tri_projected.vertex[0] = tri_projected.vertex[0] / tri_projected.vertex[0].w;
					tri_projected.vertex[1] = tri_projected.vertex[1] / tri_projected.vertex[1].w;
					tri_projected.vertex[2] = tri_projected.vertex[2] / tri_projected.vertex[2].w;

					//uninvert X/Y
					tri_projected.vertex[0].x *= -1.0f;
					tri_projected.vertex[1].x *= -1.0f;
					tri_projected.vertex[2].x *= -1.0f;
					tri_projected.vertex[0].y *= -1.0f;
					tri_projected.vertex[1].y *= -1.0f;
					tri_projected.vertex[2].y *= -1.0f;

					//offset into visible normalized space
					vec3d view_offset = { 1,1,0 };
					tri_projected.vertex[0] = tri_projected.vertex[0] + view_offset;
					tri_projected.vertex[1] = tri_projected.vertex[1] + view_offset;
					tri_projected.vertex[2] = tri_projected.vertex[2] + view_offset;
					tri_projected.vertex[0].x *= 0.5f * (float)ScreenWidth();
					tri_projected.vertex[0].y *= 0.5f * (float)ScreenHeight();
					tri_projected.vertex[1].x *= 0.5f * (float)ScreenWidth();
					tri_projected.vertex[1].y *= 0.5f * (float)ScreenHeight();
					tri_projected.vertex[2].x *= 0.5f * (float)ScreenWidth();
					tri_projected.vertex[2].y *= 0.5f * (float)ScreenHeight();

					//store poly to sort and clip against the screen edges
					polys_to_rasterize.push_back(tri_projected);
				}
			}
		}

		//painter's algorithm: take the midpoint of each triangle to sort them from back to front
		/*if (!loaded_mesh.textured)
        {
            sort(polys_to_rasterize.begin(), polys_to_rasterize.end(), [](poly &t1, poly &t2)
            {
                float z1 = (t1.vertex[0].z + t1.vertex[1].z + t1.vertex[2].z) / 3.0f;
                float z2 = (t2.vertex[0].z + t2.vertex[1].z + t2.vertex[2].z) / 3.0f;
                return z1 > z2;
            });
        }*/

		//refresh screen & clear depth buffer
		Clear(bg_color);
		for (int i = 0; i < ScreenWidth()*ScreenHeight(); i++) depth_buffer[i] = 0.0f;
		if (loaded_mesh.AIPath)
            FillRect(0, (int) ScreenHeight() * 0.5f, ScreenWidth(), (int) ScreenHeight() * 0.5f, fl_color);

		//render all transformed, viewed, projected, and sorted polys
		for (auto &tri : polys_to_rasterize)
		{
			//clip every poly against the four screen edges
			poly clipped[2];
			list<poly> poly_queue;

			poly_queue.push_back(tri);
			int num_new_polys = 1;

			//sequentially clip every triangle against each plane one at a time
			for (int step = 0; step < 4; step++)
			{
				int num_tris_to_add = 0;
				while (num_new_polys > 0)
				{
					poly test = poly_queue.front();
					poly_queue.pop_front();
					num_new_polys--;

                    //clipping based on projected viewspace rather than a viewing frustrum before projection... probably because it's easier?
					switch (step)
					{
                        case 0:	num_tris_to_add = test.ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, clipped[0], clipped[1]); break; 						//top
                        case 1:	num_tris_to_add = test.ClipAgainstPlane({ 0.0f, (float)ScreenHeight() - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, clipped[0], clipped[1]); break; 	//bottom
                        case 2:	num_tris_to_add = test.ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, clipped[0], clipped[1]); break; 						//left
                        case 3:	num_tris_to_add = test.ClipAgainstPlane({ (float)ScreenWidth() - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, clipped[0], clipped[1]); break; 	//right
					}

					//add any new triangles to the back of the queue to test against the next plane(s)
					for (int t = 0; t < num_tris_to_add; t++)
						poly_queue.push_back(clipped[t]);
				}
				num_new_polys = poly_queue.size();
			}

			//draw all the polys after clipping
			for (auto &t : poly_queue)
			{
				//if (!loaded_mesh.textured)
				//	FillTriangle(t.vertex[0].x, t.vertex[0].y, t.vertex[1].x, t.vertex[1].y, t.vertex[2].x, t.vertex[2].y, t.color);
				//else
					DrawTextured(t);
			}
		}
		return true;
	}

    void DrawTextured(poly p)
	{
		//convert 3d vertices to integers to draw on a pixel-by-pixel basis
		int x1 = (int) p.vertex[0].x, y1 = (int) p.vertex[0].y;
		int x2 = (int) p.vertex[1].x, y2 = (int) p.vertex[1].y;
		int x3 = (int) p.vertex[2].x, y3 = (int) p.vertex[2].y;
		float u1 = p.texture_vertex[0].u, v1 = p.texture_vertex[0].v, w1 = p.texture_vertex[0].w;
		float u2 = p.texture_vertex[1].u, v2 = p.texture_vertex[1].v, w2 = p.texture_vertex[1].w;
		float u3 = p.texture_vertex[2].u, v3 = p.texture_vertex[2].v, w3 = p.texture_vertex[2].w;

		//sort variables based on which vertex is higher than the others
		if (y2 < y1)
			{swap(y1, y2); swap(x1, x2); swap(u1, u2); swap(v1, v2); swap(w1, w2);}
		if (y3 < y1)
			{swap(y1, y3); swap(x1, x3); swap(u1, u3); swap(v1, v3); swap(w1, w3);}
		if (y3 < y2)
			{swap(y2, y3); swap(x2, x3); swap(u2, u3); swap(v2, v3); swap(w2, w3);}

		//y1 is the highest point, y2 is the mid point, and y3 is the lowest point
		//therefore if the triangle lays horizontally flat, dy3 is 0 and the bottom "half" does not need to be calculated
		//conversely if there is a horizontal side on the top, dy1 is 0 and only the bottom "half" needs to be calculated
		//dy2 should basically never be 0 because that means the entire triangle is just a horizontal line...
		//in all other cases, the triangle will be split into top and bottom halves
		int dy1 = y2-y1; int dx1 = x2-x1; float du1 = u2-u1; float dv1 = v2-v1; float dw1 = w2-w1;
		int dy2 = y3-y1; int dx2 = x3-x1; float du2 = u3-u1; float dv2 = v3-v1; float dw2 = w3-w1;
		int dy3 = y3-y2; int dx3 = x3-x2; float du3 = u3-u2; float dv3 = v3-v2; float dw3 = w3-w2;

		float tex_u, tex_v, tex_w;

		float dax_step = 0, dbx_step = 0,
			  du1_step = 0, dv1_step = 0, dw1_step = 0,
			  du2_step = 0, dv2_step = 0, dw2_step = 0,
			  du3_step = 0, dv3_step = 0, dw3_step = 0;

		//set steps only if the corresponding half exists
		if (dy1)
			{float fabsdy1 = (float)abs(dy1); dax_step = dx1/fabsdy1; du1_step = du1/fabsdy1; dv1_step = dv1/fabsdy1; dw1_step = dw1/fabsdy1;}
		if (dy2)
			{float fabsdy2 = (float)abs(dy2); dbx_step = dx2/fabsdy2; du2_step = du2/fabsdy2; dv2_step = dv2/fabsdy2; dw2_step = dw2/fabsdy2;}
		if (dy3)
			{float fabsdy3 = (float)abs(dy3); du3_step = du3/fabsdy3; dv3_step = dv3/fabsdy3; dw3_step = dw3/fabsdy3;}

		//scanline fill triangles pixel by pixel
		//start with the half (or entire triangle) with the flat bottom
		if (dy1)
		{
			//draw rows from top to bottom
			for (int i = y1; i <= y2; i++)
			{
				//set starting vertices
				int ax = x1 + (float)(i - y1) * dax_step;
				int bx = x1 + (float)(i - y1) * dbx_step;

				float tex_u_start = u1 + (float)(i - y1) * du1_step;
				float tex_v_start = v1 + (float)(i - y1) * dv1_step;
				float tex_w_start = w1 + (float)(i - y1) * dw1_step;

				float tex_u_end = u1 + (float)(i - y1) * du2_step;
				float tex_v_end = v1 + (float)(i - y1) * dv2_step;
				float tex_w_end = w1 + (float)(i - y1) * dw2_step;

				//make sure pixels are drawn from left to right
				if (ax > bx)
				{
					swap(ax, bx);
					swap(tex_u_start, tex_u_end);
					swap(tex_v_start, tex_v_end);
					swap(tex_w_start, tex_w_end);
				}

				tex_u = tex_u_start;
				tex_v = tex_v_start;
				tex_w = tex_w_start;

				float t_step = 1.0f / ((float)(bx - ax));
				float t = 0.0f;

				for (int j = ax; j < bx; j++)
				{
					tex_u = (1.0f - t) * tex_u_start + t * tex_u_end;
					tex_v = (1.0f - t) * tex_v_start + t * tex_v_end;
					tex_w = (1.0f - t) * tex_w_start + t * tex_w_end;
					if (tex_w > depth_buffer[i*ScreenWidth() + j])
					{
					    if (loaded_mesh.textured)
                        {
                            olc::Pixel sample = texture_spr.Sample(tex_u / tex_w, tex_v / tex_w);
                            sample *=  p.luminance;
                            Draw(j, i, sample);
                        }
                        else
                            Draw(j, i, p.color*p.luminance);
						depth_buffer[i*ScreenWidth() + j] = tex_w;
					}
					t += t_step;
				}

			}
		}

		if (dy3)
			{float fabsdy3 = (float)abs(dy3); dax_step = dx3/fabsdy3;}

		if (dy3)
		{
			for (int i = y2; i <= y3; i++)
			{
				int ax = x2 + (float)(i - y2) * dax_step;
				int bx = x1 + (float)(i - y1) * dbx_step;

				float tex_u_start = u2 + (float)(i - y2) * du3_step;
				float tex_v_start = v2 + (float)(i - y2) * dv3_step;
				float tex_w_start = w2 + (float)(i - y2) * dw3_step;

				float tex_u_end = u1 + (float)(i - y1) * du2_step;
				float tex_v_end = v1 + (float)(i - y1) * dv2_step;
				float tex_w_end = w1 + (float)(i - y1) * dw2_step;

				if (ax > bx)
				{
					swap(ax, bx);
					swap(tex_u_start, tex_u_end);
					swap(tex_v_start, tex_v_end);
					swap(tex_w_start, tex_w_end);
				}

				tex_u = tex_u_start;
				tex_v = tex_v_start;
				tex_w = tex_w_start;

				float t_step = 1.0f / ((float)(bx - ax));
				float t = 0.0f;

				for (int j = ax; j < bx; j++)
				{
					tex_u = (1.0f - t) * tex_u_start + t * tex_u_end;
					tex_v = (1.0f - t) * tex_v_start + t * tex_v_end;
					tex_w = (1.0f - t) * tex_w_start + t * tex_w_end;

					if (tex_w > depth_buffer[i*ScreenWidth() + j])
					{
						if (loaded_mesh.textured)
                        {
                            olc::Pixel sample = texture_spr.Sample(tex_u / tex_w, tex_v / tex_w);
                            sample *=  p.luminance;
                            Draw(j, i, sample);
                        }
                        else
                            Draw(j, i, p.color*p.luminance);
						depth_buffer[i*ScreenWidth() + j] = tex_w;
					}
					t += t_step;
				}
			}
		}
	}

};
#endif

