#define pi 3.1415926
struct Sphere {
    vec4 color;
    vec3 center;
    float radius;
    float reflection;
    float refractivity;
};

struct Camera {
    vec3 position;
    vec3 direction;

};

struct Ray {
    vec3 origin;
    vec3 direction;
};

struct HitRecord {
    bool hit;
    float t;
    float reflection;
    vec3 position;
    vec3 normal;
    vec3 color;
    int target;//0 is triangle;1 is sphere;2 is light
    float refractivity;
};

struct Vertex{
    vec3 position;
    vec3 normal;
    vec3 color;
};

struct Triangle{
    Vertex vertices[3];
    float reflection;
};


HitRecord intersectSphere(Ray ray, Sphere sphere) {
    HitRecord record;
    record.hit = false;
    record.t = -1.0;
    record.reflection = sphere.reflection;
    record.position = vec3(0.0);
    record.normal = vec3(0.0);
    record.color = vec3(0.0);
    record.target = -1;
    record.refractivity=-1;

    vec3 oc = ray.origin - sphere.center;
    float a = dot(ray.direction, ray.direction);
    float b = 2.0 * dot(oc, ray.direction);
    float c = dot(oc, oc) - sphere.radius * sphere.radius;
    float discriminant = b * b - 4.0 * a * c;

    if (discriminant > 0.0) {
        float t1 = (-b - sqrt(discriminant)) / (2.0 * a);
        float t2 = (-b + sqrt(discriminant)) / (2.0 * a);

        if (t1 > 0.0 || t2 > 0.0) {
            float t = min(t1, t2);
            record.hit = true;
            record.t = t;
            record.position = ray.origin + t * ray.direction;
            record.normal = normalize(record.position - sphere.center);
            record.color = sphere.color.xyz;
            record.target = 1;
            record.refractivity=sphere.refractivity;
        }
    }

    return record;
}



vec3 getBary(Triangle tri, vec3 position){

    vec3 a=tri.vertices[0].position;
    vec3 b=tri.vertices[1].position;
    vec3 c=tri.vertices[2].position;

    float detA = dot(cross(a, b), c);
    
    vec3 result;
    
    float detX = dot(cross(position, b), c);
    float detY = dot(cross(a, position), c);
    float detZ = dot(cross(a, b), position);
    
    result.x = detX / detA;
    result.y = detY / detA;
    result.z = detZ / detA;
    
    return result;
}

HitRecord intersectTriangle(const Ray ray, const Triangle tri) {
    HitRecord hitRecord;
    hitRecord.hit = false;
    hitRecord.t=-1;
    hitRecord.reflection = -1;
    hitRecord.position=vec3(-1);
    hitRecord.normal=vec3(0);
    hitRecord.color = vec3(1);
    hitRecord.target = -1;

    const vec3 e1 = tri.vertices[1].position - tri.vertices[0].position;
    const vec3 e2 = tri.vertices[2].position - tri.vertices[0].position;
    const vec3 h = cross(ray.direction, e2);
    const float a = dot(e1, h);

    if (a > -0.00001 && a < 0.00001) {
        return hitRecord;
    }

    const float f = 1.0 / a;
    const vec3 s = ray.origin - tri.vertices[0].position;
    const float u = f * dot(s, h);

    if (u < 0.0 || u > 1.0) {
        return hitRecord;
    }

    const vec3 q = cross(s, e1);
    const float v = f * dot(ray.direction, q);

    if (v < 0.0 || u + v > 1.0) {
        return hitRecord;
    }

    const float t = f * dot(e2, q);

    if (t > 0.00001) {
        hitRecord.hit = true;
        hitRecord.t=t;
        hitRecord.reflection=tri.reflection;
        hitRecord.position = ray.origin + ray.direction * t;

        vec3 bary=getBary(tri, hitRecord.position);

        hitRecord.normal = normalize(bary.x*tri.vertices[0].normal+bary.y*tri.vertices[1].normal+bary.z*tri.vertices[2].normal);
        hitRecord.color = (bary.x*tri.vertices[0].color+bary.y*tri.vertices[1].color+bary.z*tri.vertices[2].color)/(bary.x+bary.y+bary.z);
        hitRecord.target = 0;
        return hitRecord;
    }

    return hitRecord;
}

int a=747796405;
int b=277803737;

float random(int seed) {
    seed  = seed * a + 1;
    int word = ((seed >> ((seed >> 28) + 4)) ^ seed) * b;
    word      = (word >> 22) ^ word;

    a+=seed/3;
    b+=seed/2;
    return float(word) / 4294967295.0f;
}

float random(int seed,float min, float max) {
    return min + (max-min)*random(seed);
}

vec3 random_point_on_sphere(vec3 center, float radius,int seed) {
    float theta = 2.0 * 3.14159265359 * random(seed,0,1);  // Azimuthal angle
    float phi = acos(2.0 * random(seed,0,1) - 1.0);        // Polar angle

    float x = radius * sin(phi) * cos(theta);
    float y = radius * sin(phi) * sin(theta);
    float z = radius * cos(phi);

    return center + vec3(x, y, z);
}

mat3 rotate(float angle, vec3 axis) {
    axis = normalize(axis);
    float c = cos(angle);
    float s = sin(angle);
    float t = 1.0 - c;

    float x = axis.x;
    float y = axis.y;
    float z = axis.z;

    return mat3(
        t * x * x + c,     t * x * y + s * z, t * x * z - s * y,
        t * x * y - s * z, t * y * y + c,     t * y * z + s * x,
        t * x * z + s * y, t * y * z - s * x, t * z * z + c
    );
}

vec3 getVectorFromAngles(float theta_x, float theta_y) {
    float x = sin(theta_y) * cos(theta_x);
    float y = cos(theta_y);
    float z = sin(theta_y) * sin(theta_x);
    return vec3(x, y, z);
}

vec3 random_reflection_direction_0(vec3 normal, vec3 lightIn, float theta,int seed) {

    float cos_phi=dot(normalize(lightIn),normalize(normal));

    float sin_theta=sin(theta);
    vec3 random_point = vec3(0);

    float random_x_theta=random(seed/2,-pi,pi);

    float random_y_theta=random(seed/2,-pi,pi);

    random_point=normalize(getVectorFromAngles(random_x_theta,random_y_theta));

    vec3 reflected_specular = normalize(reflect(lightIn, normal));
    vec3 result=reflected_specular+sin_theta*random_point;
    if(dot(result,normal)<0){
        return vec3(0);
    }
    return result;
    
}
vec3 random_reflection_direction_1(vec3 normal, vec3 lightIn, float theta,int seed) {

    float random_x_theta=random(seed/2,-pi,pi);

    float random_y_theta=random(seed/2,-pi,pi);


    vec3 random=getVectorFromAngles(random_x_theta,random_y_theta);


    

    return normalize(normal)+random*sin(theta);

}


mat3 compute_transformation_matrix(vec3 a, vec3 b) {
    vec3 c = normalize(cross(a, b));
    vec3 d = normalize(cross(c, a));
    return mat3(a, c, d);
}


vec3 compute_refracted_ray(vec3 lightIn, vec3 normal, float n) {
    // 计算入射角的余弦值
    float cosTheta1 = dot(lightIn, normal);

    // 使用斯涅尔定律计算折射角的正弦值
    float sinTheta2 = sin(acos(cosTheta1)) / n;

    // 计算折射光线的方向
    vec3 refractedDir = normalize(n * lightIn - (n * cosTheta1 + sqrt(1.0 - n * n * (1.0 - cosTheta1 * cosTheta1))) * normal);

    return refractedDir;
}

vec3 perpendicular_unit_vector(vec3 p, vec3 target) {

    vec3 projection = dot(p, target) * target;

    vec3 perpendicular = p - projection;

    return normalize(perpendicular);
}

vec3 symmetric_point(vec3 p, vec3 axis) {
    vec3 projection = dot(p, axis) * axis;
    vec3 displacement = 2.0 * projection;
    vec3 symmetricPoint = p - displacement;
    return symmetricPoint;
}

bool is_total_internal_reflection(vec3 lightIn, vec3 normal, float n) {
    // 计算入射角的余弦值
    float cosTheta1 = dot(lightIn, normal);
    
    // 计算临界角的正弦值
    float sinCriticalAngle = n;
    
    // 判断是否发生全反射
    return cosTheta1 > sinCriticalAngle;
}



vec3 get_refract_direction__in_a_sphere(vec3 normal, vec3 lightIn, float n,vec3 center){

    if(is_total_internal_reflection(lightIn,normal,n)){
        return reflect(lightIn,normal);
    }

    vec3 refrection=compute_refracted_ray(lightIn,normal,n);
    vec3 like_a_normal=perpendicular_unit_vector(center,refrection);

    vec3 lightOut=reflect(lightIn,like_a_normal);

    if(dot(lightIn,lightOut)<0){
        lightOut=-lightOut;
    }
    return lightOut;

}

vec3 get_refract_position__in_a_sphere(vec3 normal, vec3 lightIn, float n,vec3 center,vec3 inPosition){

    if(is_total_internal_reflection(lightIn,normal,n)){
        return inPosition;
    }

    vec3 refrection=compute_refracted_ray(lightIn,normal,n);
    vec3 like_a_normal=perpendicular_unit_vector(center,refrection);

    return symmetric_point(inPosition,like_a_normal);
}

vec3 compute_vector_c(vec3 a, vec3 b, float dot_bc) {
    float k1 = dot_bc / dot(b, a);
    float k2 = (dot_bc - k1 * dot(b, a)) / dot(b, b);
    
    vec3 vector_c = k1 * a + k2 * b;
    return vector_c;
}

