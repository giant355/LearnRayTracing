#include "geometry.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include "tgaimage.h"
#include <iostream>
#include <vector>

const int width = 800;
const int height = 600;
const double fov = 3.1415926 / 3.5;

struct Material {
    TGAColor color;
};

struct Light {
    Light(const vec3& p, const float& i) : position(p), intensity(i) {}
    vec3 position;
    float intensity;
};

struct Sphere {
    vec3 center;
    double r;
    Material material;

    double ray_intersect(const vec3& orig, const vec3& dir,vec3& hit_point) const {
        vec3 D = normalized(dir);
        vec3 L = orig - center;
        double a = 1,
               b = 2 * (L * D),
               c = L * L - r * r;
        double delta = b * b - 4 * a * c;
        if (delta >= 0) {
            double t1 = (-b - sqrt(delta)) / (2 * a);
            double t2 = (-b + sqrt(delta)) / (2 * a);
            
            if (t1 < 0)std::swap(t1, t2);
            if (t1 < 0)return -1.;

            hit_point = orig + t1 * D;//此时t1>=0;
            return t1;
        }
        return -1;//delta<0
    }
}; 

bool scene_intersect(const vec3& orig, const vec3& dir, const std::vector<Sphere>& spheres, const std::vector<Light>& lights,vec3& hit_point,vec3& N,Material& material) {
    double t = std::numeric_limits<double>::max();
    for (const Sphere& sphere : spheres) {
        vec3 hitPos;
        double t1 = sphere.ray_intersect(orig, dir, hitPos);
        if (t1 < 0)continue;
        if (t1 < t) {
            t = t1;
            hit_point = hitPos;
            N = normalized(hit_point - sphere.center);
            material.color = sphere.material.color;
        }
    }
    return t < 1000.;
}

vec3 reflect(const vec3& I, const vec3& N) {
    return (I - 2 * (N * I) * N);
}

TGAColor ray_cast(const vec3& orig,const vec3& dir,const std::vector<Sphere>& spheres,const std::vector<Light>& lights,int depth = 0) {
    vec3 hit_point, N;
    Material material;
    double diffuse = 0., specular = 0.;
    TGAColor reflect_color;
   
    if (!scene_intersect(orig,dir,spheres,lights,hit_point,N,material)||depth >= 4)return { 10,10,10 };//啥也每打到，返回背景色

    vec3 reflectDir = reflect(dir, N);
    vec3 reflectOrig = (reflectDir * N < 0) ? (hit_point - reflectDir * 1e-3) : (hit_point + reflectDir * 1e-3);
    reflect_color = ray_cast(reflectOrig, reflectDir, spheres, lights, depth + 1);

    for (const Light& light : lights) {
        vec3 lightDir = normalized(light.position - hit_point);
        vec3 shadowOrig = (lightDir * N < 0) ? (hit_point - N * 1e-3) : (hit_point + N * 1e-3);
        double lightDistance = norm(light.position - shadowOrig);
        vec3 shadow_pt, shadow_N; Material trash;

        if (scene_intersect(shadowOrig, lightDir, spheres, lights, shadow_pt, shadow_N, trash)&& norm(shadow_pt - shadowOrig) < lightDistance) continue;

        diffuse += std::max(0., lightDir * N) * light.intensity;
        vec3 r = normalized(2 * (N * lightDir) * N - lightDir);
        vec3 eye = normalized(orig - hit_point);
        specular += std::pow(std::max(eye * r, 0.), 35) * light.intensity;
    }
    for (int channel : {0, 1, 2}) {
        material.color[channel] = std::min(255.,material.color[channel] * (diffuse*.4+specular*2.)+reflect_color[channel]*1.5);
    }
    return material.color;
}

void rend(TGAImage& framebuffer, const std::vector<Sphere>& spheres,const std::vector<Light>& lights) {
    vec3 orig{ 0, 0, 0 };

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            double x = (2 * (i + 0.5) / static_cast<double>(width) - 1) * tan(fov / 2) * width / static_cast<double>(height);
            double y = (2 * (j + 0.5) / static_cast<double>(height) - 1) * tan(fov / 2);

            framebuffer.set(i, j, ray_cast(orig, {x,y,-1},spheres,lights));
        }
    }

    framebuffer.write_tga_file("framebuffer.tga");
}

int main(int argc, char** argv) {
    TGAImage framebuffer(width, height, TGAImage::RGB);

    Material ivory       {{102, 102, 77, 255}};
    Material red_rubber  {{77, 26, 26, 255}};

    std::vector<Sphere> objs{
        {{-3, 0, -16}, 2., ivory},
        {{-1.0, -1.5, -12}, 2., ivory},
        {{1.5, -0.5, -18}, 3., red_rubber},
        {{7, 5, -18}, 4., ivory},
    };

    std::vector<Light> lights{
       {{-20, 20, 20}, 1.5},
       {{30, 50, -25}, 3.},
       {{30, 20, 30}, 2.},
    };

    rend(framebuffer, objs,lights);
    std::cout << "Render complete: framebuffer.tga" << std::endl;
    return 0;
}
