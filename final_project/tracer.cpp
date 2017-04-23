Vec3f trace(
    const Vec3f &rayorigin, 
    const Vec3f &raydir, 
    const std::vector<Sphere> &spheres, 
    const int &depth)
{
    float tnear = INFINITY; 
    const Sphere* sphere = NULL; // This sphere is for storing the intersected spphere
    // find intersection of this ray with the spehre in the scene 
    // Traversing through all the spheres. 
    for(unsigned i = 0; i<spheres.size(); i++){
        float t0 = INFINITY, t1 = INFINITY; 
        if(spheres[i].intersect(rayorig, raydir, t0, t1)){
            if (t0 < 0) t0= t1; 
            if (t0 < tnear){
                tnear = t0; 
                sphere = &spheres[i]; 
            }
        }
    }

    // If there's no intersection return black or background color 
    if (!sphere) return Vec3f(2); 
    Vec3f surfaceColor = 0; 
    Vec3f phit = rayorig + raydir * tnear; // point of intersection 
    Vec3f nhit = phit - sphere->center;    // normal at the intersection point 
    nhit.normalize();                      // normalize normal direction 

    if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
    if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH){
        float facingratio = -raydir.dot(nhit); 

        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1); 

        // compute reflection direction 
        Vec3f refldir = raydir - nhit*2*raydir.dot(nhit); 
        refldir.normalize(); 
        Vec3f reflection = trans(phit + nhit * bias, refldir, spheres, depth + 1); 
        Vec3f refraction = 0; 

        // If the sphere is also transparent compute refraction ray (transmission)
        if(sphere->transparency){
            float ior 1.1, eta = (inside) ? ior : 1/ior; // are we inside or outside the surface
            float cosi ...... 
        }
    }
}

/***** SUDO CODE ****/ 
// A function that given the rayorigin and the ray direction, computes the first 
// hit position (new origin) and then the refraction direction. 
// This should be the function that gets called recursively. 

// given number of slabs 

// while this function does not intersect with a sphere, then keep calling it. 
bool nextIntersection(
    const Vec3f &rayorigin, 
    const Vec3f &raydir, 
    Vec3f &new_origin, 
    Vec3f &new_direction, 
    const std::vector<Sphere> &spheres, 
    const int &depth, 
    int num_slabs)
{
    int num_slabs = 5;   
    double height = 480; // NEED TO CHANGE THIS VALUE
    double slab_height = height / num_slabs;

    // check if there is intersection 
    for(unsigned i = 0; i<spheres.size(); i++){
        if (spheres[i].intersect(rayorigin, raydir, t0, t1)) return false; 
    }

    double alpha = 0;
    double ceiling, floor;  
    if (raydir.y >= 0){
        alpha = floor(rayorigin.y / slab_height + slab_height);
        ceiling = alpha; 
        floor = alpha - slab_height;  
    } 
    else{
        alpha = floor(rayorigin.y / slab_height);
        ceiling = alpha + slab_height; 
        floor = alpha;  
    }

 	double index = linear_refraction_interpolate((ceiling + floor)/2.0, 1, 1.05); 
    double t = (alpha - rayorigin.x)/(raydir.y); 
    
    new_origin.x = rayorigin.x + t*raydir.x; 
    new_origin.y = rayorigin.y + t*raydir.y; 
    new_origin.z = rayorigin.z + t*raydir.z; 

    // calculating the new direction 
    double theta_1 = 180 - acos(dot(raydir, Vector3D(0, 1, 0)));
    double theta_2 = asin(index1/index2 * sin(theta_1));  
    new_direction = Vector3D(0, -1, 0)/ cos(theta_2); 

    return true; 
}