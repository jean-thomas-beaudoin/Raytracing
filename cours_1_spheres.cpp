#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#include <iostream>
#include <limits>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.1415926535897

// définition de la classe vecteur pour les vecteurs et les points

class vector{
	public :
		vector(double x = 0, double y = 0, double z = 0){
            coord[0] = x;
            coord[1] = y;
            coord[2] = z;
        };
        
        const double& operator[](int i) const {return coord[i];}
        double operator[](int i){return coord[i];}  // récupération des éléments des vecteurs via leur position (0, 1 ou 2)

        double norm2(){
            return coord[0]*coord[0]+coord[1]*coord[1]+coord[2]*coord[2];  // norme du vecteur
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

// surcharge des opérateurs pour faire les opérations avec les vecteurs

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

// définition de la classe des rayons

class Ray{
	public:
		Ray(const vector& C, const vector& U) : c(C), u(U) {};
		vector c,u ;
};

// définition de la classe des sphères

class sphere{
	public:
		sphere(const vector& O, double R, const vector& color) : O(O), R(R), albedo(color) {};
        vector O;
        double R;
        vector albedo;
        bool intersect(Ray& d, vector& pos, vector& n, double& t){
            double b = 2*dot(d.u,d.c-O);
            double c = (d.c-O).norm2()-R*R;
            double delta = b*b-4*c;
            if(delta < 0){
                return false;
            };
            double t1 = (-b-sqrt(delta))/2;
            double t2 = (-b+sqrt(delta))/2;
            if (t2 < 0){
                return false;
            } else {
                if (t1 > 0){
                    t = t1;
                } else {
                    t = t2;
                };
                pos = d.c + t*d.u;
                n = (pos - O).getnormalized();

                return true;
            };  
        };
};

// définition de la classe scene contenant toutes les spheres étant définies pour la réalisation de l'image

class scene{
    public :
        scene(){};
        void add_sphere(const sphere& s){
            spheres.push_back(s);
        };
    std::vector<sphere> spheres;

        // fonction intersect qui renvoie true si il y a intersection entre un rayon et une sphere de la scene et false sinon

        bool intersect(Ray& d, vector& pos, vector& n, int& sphere_id){

            bool has_intersection = false;
            double t_min = std::numeric_limits<double>::max();
            sphere_id = -1;
            for (int i = 0 ; i < spheres.size(); i++){
                vector local_pos, local_n;
                double t;
                bool has_intersection_local = spheres[i].intersect(d, local_pos, local_n, t);
                if (has_intersection_local){
                    has_intersection = true;
                    if (t < t_min){
                        t_min = t;
                        pos = local_pos;
                        n = local_n;
                        sphere_id = i;
                    };
                };
            };
        if (has_intersection && sphere_id == -1){
            std::cout << "bug" << std::endl;
        };
        return has_intersection;
        };

};

// fonctions max et min

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

// programme général pour la réalisation de l'image

int main() {
	int H = 512;    
    int W = 512;
    double fov = 60*M_PI/180;
	
    sphere s1(vector (0,0,0) ,10 , vector (1,1,1)); // sphere centrale blanche
    sphere s2(vector (0,-1000,0), 990, vector(0,0,1)); // sol bleu
    sphere s3(vector (0,1000,0), 990, vector(1,0,0)); // plafond rouge
    //sphere s4(vector (0,0,-1000), 990, vector(0,1,0)); // droite verte
    //sphere s5(vector (0,0,1000), 990, vector(1,0,1)); // gauche violette

    scene s;
    s.add_sphere(s1);
    s.add_sphere(s2);
    s.add_sphere(s3);
    //s.add_sphere(s4);
    //s.add_sphere(s5);

    vector C(0.,0.,55.); // position de la caméra (point de départ du rayon)

    vector pos_light(15,60,40); // position de la lumière
    double intensity_light = 1000000;

	std::vector<unsigned char> image(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
            vector u(j-W/2,i-H/2,-W/(2*tan(fov/2))); // direction du rayon partant de C
            u.normalize();

            Ray r(C,u);
            vector pos, n;
            int sphere_id;
            bool has_intersection = s.intersect(r,pos,n,sphere_id);
            double intensity_pix = 0;
            if (has_intersection){

                intensity_pix = intensity_light * max(dot((pos_light - pos).getnormalized(),n) ,0)/ (pos_light - pos).norm2();
                vector intensity_albedo = intensity_pix * s.spheres[sphere_id].albedo;
			    image[((H-i-1)*W + j) * 3 + 0] = min(255, max(intensity_albedo[0],0));
			    image[((H-i-1)*W + j) * 3 + 1] = min(255, max(intensity_albedo[1],0));
			    image[((H-i-1)*W + j) * 3 + 2] = min(255, max(intensity_albedo[2],0));
            };
		};
	};
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
    

	return 0;
};