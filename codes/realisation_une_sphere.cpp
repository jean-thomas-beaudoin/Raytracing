#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.1415926535897

class vector{
	public :
		vector(double x = 0, double y = 0, double z = 0){
            coord[0] = x;
            coord[1] = y;
            coord[2] = z;
        };
        
        const double& operator[](int i) const {return coord[i];}
        double operator[](int i){return coord[i];}

        double norm2(){
            return coord[0]*coord[0]+coord[1]*coord[1]+coord[2]*coord[2];
        };

        void normalize(){
            double n = sqrt(norm2());
            coord[0] /= n;
            coord[1] /= n;
            coord[2] /= n;
        };

        vector getnormalized(){
            vector result(*this);
            result.normalize();
            return result;
        };

	private :
		double coord[3];
};

vector operator+(const vector& a, const vector& b){
    return vector(a[0]+b[0],a[1]+b[1],a[2]+b[2]);
};

vector operator-(const vector& a, const vector& b){
    return vector(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
};

vector operator*(double a, const vector& b){
    return vector(a*b[0],a*b[1],a*b[2]);
};

vector operator*(const vector& b, double a){
    return vector(a*b[0],a*b[1],a*b[2]);
};

vector operator*(const vector& a, const vector& b){
    return vector(a[0]*b[0],a[1]*b[1],a[2]*b[2]);
};

vector operator/(const vector& b, double a){
    return vector(b[0]/a,b[1]/a,b[2]/a);
};

double dot(const vector& a, const vector& b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
};

class Ray{
	public:
		Ray(const vector& C, const vector& U) : c(C), u(U) {};
		vector c,u ;
};

class sphere{
	public:
		sphere(const vector& O, double R, const vector& color) : O(O), R(R), albedo(color) {};
        vector O;
        double R;
        vector albedo;
};



double max(double a, double b){
    if (a > b){
        return a;
    } else {
        return b;
    };
};

double min(double a, double b){
    if (a < b){
        return a;
    } else {
        return b;
    };
};

bool intersect(Ray& d, sphere& s, vector& pos, vector& n){
    double b = 2*dot(d.u,d.c-s.O);
    double c = (d.c-s.O).norm2()-s.R*s.R;
    double delta = b*b-4*c;
    if(delta < 0){
        return false;
    };
    double t1 = (-b-sqrt(delta))/2;
    double t2 = (-b+sqrt(delta))/2;
    if (t2 < 0){
        return false;
    } else {
        double t;
        if (t1 > 0){
            t = t1;
        } else {
            t = t2;
        };
        pos = d.c + t*d.u;
        n = (pos - s.O).getnormalized();

        return true;
    };  
};


int main() {
	int H = 512;    
    int W = 512;
    double fov = 60*M_PI/180;
	
    vector O(0.,0.,0.);
    int R = 10;
    sphere s(O,R, vector (1,1,0));
    vector C(0,0,55);

    vector pos_light(15, 100,40);
    double intensity_light = 6000000;

	std::vector<unsigned char> image(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
            vector u(j-W/2,i-H/2,-W/(2*tan(fov/2)));
            u.normalize();

            Ray r(C,u);
            vector pos, n;
            bool has_intersection = intersect(r,s,pos,n);
            double intensity_pix = 0;
            if (has_intersection){
                intensity_pix = intensity_light * max(dot((pos_light - pos).getnormalized(),n) ,0)/ (pos_light - pos).norm2();
                vector intensity_albedo = s.albedo * intensity_pix;
			    image[((H-i-1)*W + j) * 3 + 0] = min(255, max(intensity_albedo[0],0));
			    image[((H-i-1)*W + j) * 3 + 1] = min(255, max(intensity_albedo[1],0));
			    image[((H-i-1)*W + j) * 3 + 2] = min(255, max(intensity_albedo[2],0));
            };
		};
	};
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
    

	return 0;
};