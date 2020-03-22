#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#include <iostream>
#include <limits>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.1415926535897

#include <random>

static std::default_random_engine engine;
static std::uniform_real_distribution<double> uniform(0,1);

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

        vector& operator+=(const vector& b){
            coord[0] += b[0];
            coord[1] += b[1];
            coord[2] += b[2];
            return *this;
}

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

vector operator-(const vector& a){
    return vector(-a[0],-a[1],-a[2]);
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

vector cross(const vector& a, const vector& b){
    return vector(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
};

vector random_cos(vector &N){
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    vector dir_aleat_repere_local(cos(2*M_PI*r1)*sqrt(1-r2),sin(2*M_PI*r1)*sqrt(1-r2),sqrt(r2));
    vector aleat(uniform(engine)-0.5,uniform(engine)-0.5,uniform(engine)-0.5);
    vector tangent1 = cross(N,aleat);
    tangent1.normalize();
    vector tangent2 = cross(tangent1,N);
    tangent2.normalize();

    return dir_aleat_repere_local[2]*N + dir_aleat_repere_local[0]*tangent1 + dir_aleat_repere_local[1]*tangent2;
};


// définition de la classe des rayons

class Ray{
	public:
		Ray(const vector& C, const vector& U) : c(C), u(U) {};
        vector c,u;
};

// définition de la classe des sphères

class object {
    public:
        object(){};
        virtual bool intersect(const Ray& d, vector& pos, vector& n, double& t) const = 0;
        vector albedo;
        bool isMirror;
        bool isTransparent;
};

class triangle : object {
public :
    triangle(const vector& a, const vector& b, const vector& c, const vector& color, const bool& miroir = false, const bool& transparent = false) : a(a), b(b), c(c) {
        albedo = color;
        isMirror = miroir;
        isTransparent = transparent;
    };
    bool intersect(const Ray& d, vector& pos, vector& n, double& t) const {
        n = cross(b-a,c-a).getnormalized();
        t = dot(c-d.c,n)/dot(d.u,n);
        if (t<0){return false;};

        pos = d.c + t*d.u;
        vector u = b-a;
        vector v = c-a;
        vector w = pos-a;
        double m11 = u.norm2();
        double m12 = dot(u,v);
        double m22 = v.norm2();
        double detm = m11*m22-m12*m12;

        double b11 = dot(w,u);
        double b21 = dot(w,v);
        double detb = b11*m22-b21*m12;
        double beta = detb/detm;

        double g12 = b11;
        double g22=b21;
        double detg = m11*g22-m12*g12;
        double gamma = detg/detm;

        double alpha = 1-beta-gamma;
        if (alpha < 0 || alpha > 1){return false;};
        if (beta < 0 || beta > 1){return false;};
        if (gamma < 0 || gamma > 1){return false;};
        if (alpha + beta + gamma > 1){return false;};
        return true;

    };
    vector a,b,c;
};

class sphere : object {
	public:
		sphere(const vector& O, double R, const vector& color, const bool& miroir = false, const bool& transparent = false) : O(O), R(R) {
            albedo = color;
            isMirror = miroir;
            isTransparent = transparent;
        };
        vector O;
        double R;
        
        bool intersect(const Ray& d, vector& pos, vector& n, double& t) const {
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

class scene {
    public :
        scene(){};
        void add_sphere(const sphere& s){
            objects.push_back((object*)&s);
        };
        void add_triangle(const triangle& s){
            objects.push_back((object*)&s);
        };
    std::vector<object*> objects;
    sphere *light;
    double intensity_light;

        // fonction intersect qui renvoie true si il y a intersection entre un rayon et une sphere de la scene et false sinon

        bool intersect(const Ray& d, vector& pos, vector& n, int& sphere_id, double& t_min) const {

            bool has_intersection = false;
            t_min = std::numeric_limits<double>::max();
            sphere_id = -1;
            for (int i = 0 ; i < objects.size(); i++){
                vector local_pos, local_n;
                double t;
                bool has_intersection_local = objects[i]->intersect(d, local_pos, local_n, t);
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


// Fonction récursive pour la programmation des miroirs

vector getColor(const Ray& r, const scene& s, int nb_reb){
    
    if (nb_reb = 0){
        return vector (0,0,0);
    };
    vector pos, n;
    int sphere_id;
    double t;
    bool has_intersection = s.intersect(r,pos,n,sphere_id,t);
    double intensity_pix = 0;
    vector intensity_albedo (0,0,0);
    if (has_intersection){

        if(s.objects[sphere_id]->isMirror){
            vector dir_mirror = r.u-2*dot(r.u,n)*n;
            Ray ray_mirror(pos + 0.001*n,dir_mirror);
            intensity_albedo = getColor(ray_mirror,s,nb_reb-1);
        } else {
            if (s.objects[sphere_id]->isTransparent){
                double n1 = 1;
                double n2 = 1.3;
                vector normal_transp = n;
                if(dot(r.u,n) > 0){
                    n1 = 1;
                    n2 = 1;
                    normal_transp = -n;
                };
                double racine = 1-(n1*n1)/(n2*n2)*(1-dot(r.u,normal_transp)*dot(r.u,normal_transp));
                if (racine < 0){
                    intensity_albedo = vector (0,0,0);
                } else {
                    vector dir_transmis = n1/n2*r.u-(n1/n2*dot(r.u,normal_transp)+sqrt(racine))*normal_transp;
                    Ray ray_transmis(pos - 0.001*normal_transp, dir_transmis);
                    intensity_albedo = getColor(ray_transmis,s,nb_reb);
                };
            } else {

                //contribution éclairage direct
                /*Ray ray_light(pos + 0.01*n,(s.pos_light-pos).getnormalized());
                vector P_light, N_light;
                int sphere_id_light;
                double t_light;
                bool has_intersection_light = s.intersect(ray_light, P_light, N_light, sphere_id_light,t_light);
                double d_light2 = (s.pos_light - pos).norm2();

                if (has_intersection_light && t_light*t_light < d_light2){
                    intensity_pix = 0;
                } else {
                intensity_pix = s.intensity_light * max(dot((s.pos_light - pos).getnormalized(),n) ,0)/ d_light2;
                };
                intensity_albedo = intensity_pix * s.spheres[sphere_id].albedo / M_PI;*/

                vector axe = (pos-s.light->O).getnormalized();
                vector dir_aleatoire = random_cos(axe);
                vector point_aleat = dir_aleatoire * s.light->R + s.light->O;
                vector wi = (point_aleat-pos).getnormalized();
                double d_light2 = (point_aleat - pos).norm2();
                vector np = dir_aleatoire;

                Ray ray_light(pos + 0.01*n,wi);
                vector P_light, N_light;
                int sphere_id_light;
                double t_light;
                bool has_intersection_light = s.intersect(ray_light, P_light, N_light, sphere_id_light,t_light);

                if (has_intersection_light && t_light*t_light < d_light2*0.99){
                    intensity_pix = 0;
                } else {
                intensity_pix = s.intensity_light / (4*M_PI*d_light2) * max(0,dot(n, wi)) * dot(np,-wi) / dot(axe,dir_aleatoire);
                };

                intensity_albedo = intensity_pix * s.objects[sphere_id]->albedo;

                //ajout de l'éclairage indirect

                vector dir_aleat = random_cos(n);
                Ray ray_aleat(pos + 0.001*n,dir_aleat);

                intensity_albedo += getColor(ray_aleat,s,nb_reb-1) * s.objects[sphere_id]->albedo;
            };
        };
    };
    return intensity_albedo;
};

// programme général pour la réalisation de l'image

int main() {
	int H = 512;    
    int W = 512;
    const int nrays = 100;
    double fov = 60*M_PI/180;

	sphere slight(vector(15,70,-30),15,vector(1,1,1));

    sphere s1(vector (-15,0,-70) ,10 , vector (1,1,1), false, true); // sphere a gauche de la figure
    sphere s1c(vector (15,0,-130) ,10 , vector (1,1,1), true); // sphere a droite de la figure
    sphere s2(vector (0,-2000-20,0), 2000, vector(0,0,0.5)); // sol bleu
    sphere s3(vector (0,2000+100,0), 2000, vector(0.28,0.24,0.55)); // plafond rouge
    sphere s4(vector (-2000-50,0,0), 2000, vector(1,0.65,0)); // droite orange
    sphere s5(vector (2000+50,0,0), 2000, vector(0.55,0,0)); // gauche violette
    sphere s6(vector (0,0,-2000-200), 2000, vector(0,0.5,0)); // fond vert

    triangle tri(vector(-10,-10,-100), vector(10,-10,-100), vector(0,10,-100), vector(1,1,1));

    scene s;
    s.add_sphere(slight);
    s.add_triangle(tri);
    s.add_sphere(s1);
    s.add_sphere(s1c);
    s.add_sphere(s2);
    s.add_sphere(s3);
    s.add_sphere(s4);
    s.add_sphere(s5);
    s.add_sphere(s6);
    s.light = &slight; // position de la lumière
    s.intensity_light = 3000000000;
    vector pos_cam(0,0,0);
    double focus_dist = 100;
    double aperture = 0.5;

    vector C(0.,0.,55.); // position de la caméra (point de départ du rayon)

	std::vector<unsigned char> image(W*H * 3, 0);

    #pragma omp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

            vector color(0,0,0);
            for(int k=0; k < nrays; k++){
                
                double r1 = uniform(engine);
                double r2 = uniform(engine);
                double R = sqrt(-2*log(r1));
                double dx = R*cos(2*M_PI*r2);
                double dy = R*sin(2*M_PI*r2);

                double dx_aperture = (uniform(engine) - 0.5) * aperture;
                double dy_aperture = (uniform(engine) - 0.5) * aperture;

                vector u(j-W/2 + 0.5 + dx ,i-H/2 + 0.5 + dy ,-W/(2*tan(fov/2))); // direction du rayon partant de C
                u.normalize();
                
                vector destination = pos_cam + focus_dist *u;
                vector new_origin = pos_cam + vector(dx_aperture,dy_aperture,0);

                Ray r(new_origin,(destination - new_origin).getnormalized());
                color += getColor(r,s,5) / nrays;
            };

            image[((H-i-1)*W + j) * 3 + 0] = min(255, max(std::pow(color[0],1/2.2),0));
		    image[((H-i-1)*W + j) * 3 + 1] = min(255, max(std::pow(color[1],1/2.2),0));
		    image[((H-i-1)*W + j) * 3 + 2] = min(255, max(std::pow(color[2],1/2.2),0));
		};
	};
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
    

	return 0;
};