// set the precision of the float values (necessary if using float)
#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif
precision mediump int;

// define constant parameters
// EPS is for the precision issue (see precept slide)
#define INFINITY 1.0e+12
#define EPS 1.0e-3
#define PI 3.14159265359

// define constants for scene setting 
#define MAX_LIGHTS 10

// define texture types
#define NONE 0
#define CHECKERBOARD 1
#define MYSPECIAL 2

// define material types
#define BASICMATERIAL 1
#define PHONGMATERIAL 2
#define LAMBERTMATERIAL 3

// define reflect types - how to bounce rays
#define NONEREFLECT 1
#define MIRRORREFLECT 2
#define GLASSREFLECT 3

struct Shape {
    int shapeType;
    vec3 v1;
    vec3 v2;
    float rad;
};

struct Material {
    int materialType;
    vec3 color;
    float shininess;
    vec3 specular;

    int materialReflectType;
    float reflectivity; 
    float refractionRatio;
    int special;

};

struct Object {
    Shape shape;
    Material material;
};

struct Light {
    vec3 position;
    vec3 color;
    float intensity;
    float attenuate;
};

struct Ray {
    vec3 origin;
    vec3 direction;
};

struct Intersection {
    vec3 position;
    vec3 normal;
};

// uniform
uniform mat4 uMVMatrix;
uniform int frame;        
uniform float height;
uniform float width;
uniform vec3 camera;
uniform int numObjects;
uniform int numLights;
uniform Light lights[MAX_LIGHTS];
uniform vec3 objectNorm;

varying vec2 v_position;

// find then position some distance along a ray
vec3 rayGetOffset( Ray ray, float dist ) {
    return ray.origin + ( dist * ray.direction );
}

// if a newly found intersection is closer than the best found so far, record the new intersection and return true;
// otherwise leave the best as it was and return false.
bool chooseCloserIntersection( float dist, inout float best_dist, inout Intersection intersect, inout Intersection best_intersect ) {
    if ( best_dist <= dist ) return false;
    best_dist = dist;
    best_intersect.position = intersect.position;
    best_intersect.normal   = intersect.normal;
    return true;
}

// put any general convenience functions you want up here
// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 135 lines of code.

// shortest distance from plane to point
float signedDistance( vec3 norm, vec3 point, vec3 p ) {
    float d = -dot(norm, point);
    return (dot(norm, p) + d) / length(norm);
}
// ----------- STUDENT CODE END ------------

// forward declaration
float rayIntersectScene( Ray ray, out Material out_mat, out Intersection out_intersect );

// Plane
// this function can be used for plane, triangle, and box
float findIntersectionWithPlane( Ray ray, vec3 norm, float dist, out Intersection intersect ) {
    float a   = dot( ray.direction, norm );
    float b   = dot( ray.origin, norm ) - dist;
    
    if ( a < 0.0 && a > 0.0 ) return INFINITY;
    
    float len = -b/a;
    if ( len < EPS ) return INFINITY;

    intersect.position = rayGetOffset( ray, len );
    intersect.normal   = norm;
    return len;
}

// find the intersection of a ray with a plane defined by 3 points
float findIntersectionWithPlane3( Ray ray, vec3 p1, vec3 p2, vec3 p3, out Intersection out_intersect) {
    vec3 norm = cross( p3 - p1, p2 - p3 );
    return findIntersectionWithPlane(ray, norm, dot(p1, norm), out_intersect);
}

// return the area of a triangle defined by 3 points
float areaOfTriangle( vec3 t1, vec3 t2, vec3 t3 ) {
    return 0.5 * length(cross(t2 - t1, t3 - t1));
}

// return true if point p is in the triangle defined by t1, t2, t3
bool pointInTriangle( vec3 p, vec3 t1, vec3 t2, vec3 t3) {

    float a = areaOfTriangle(p, t1, t2) / areaOfTriangle(t1, t2, t3);
    float b = areaOfTriangle(p, t1, t3) / areaOfTriangle(t1, t2, t3);
    float c = areaOfTriangle(p, t2, t3) / areaOfTriangle(t1, t2, t3);

    if (a >= -EPS && a <= 1.0+EPS && b >= -EPS && b <= 1.0+EPS && c >= -EPS && c <= 1.0+EPS) return true;
    else return false;
}

// Triangle
float findIntersectionWithTriangle( Ray ray, vec3 t1, vec3 t2, vec3 t3, out Intersection intersect ) {

    // find intersection of ray with triangle plane
    float len = findIntersectionWithPlane3(ray, t1, t2, t3, intersect);
    vec3 p = intersect.position;

    // check if point is inside triangle algebraically
    vec3 v1 = t1 - ray.origin;
    vec3 v2 = t2 - ray.origin;
    vec3 v3 = t3 - ray.origin;
    vec3 n1 = normalize(cross(v2, v1));
    vec3 n2 = normalize(cross(v3, v2));
    vec3 n3 = normalize(cross(v1, v3));
    if (signedDistance(n1, ray.origin, p) < EPS || signedDistance(n2, ray.origin, p) < EPS
        || signedDistance(n3, ray.origin, p) < EPS) return INFINITY;

    v1 = t1 - p;
    v2 = t2 - p;
    v3 = t3 - p;
    n1 = normalize(cross(v2, v1));
    n2 = normalize(cross(v3, v2));
    n3 = normalize(cross(v1, v3));
    if (dot(ray.direction, n1) < EPS || dot(ray.direction, n2) < EPS || dot(ray.direction, n3) < EPS)
        return INFINITY;
   
    return len;
}

// Sphere
float findIntersectionWithSphere( Ray ray, vec3 center, float radius, out Intersection intersect ) {   
  
    vec3 L = center - ray.origin;
    float t_ca = dot( L, ray.direction );
    if (t_ca < -EPS) return INFINITY;
    float d = dot(L, L) - t_ca * t_ca;
    if (d > radius*radius + EPS) return INFINITY;
    float t_hc = sqrt(radius*radius - d);
    float t1 = t_ca - t_hc;
    float t2 = t_ca + t_hc;

    float t;

    if (t1 > EPS) t = t1;
    else if (t2 > EPS) t = t2;
    else return INFINITY;

    intersect.position = rayGetOffset(ray, t);
    intersect.normal = normalize(intersect.position - center);
    return length(t * ray.direction);
}

// Box
float findIntersectionWithBox( Ray ray, vec3 pmin, vec3 pmax, out Intersection out_intersect ) {

    float len;
    Intersection pClosest;
    float minDist = INFINITY;
    vec3 norm;
    vec3 p;

    // check each face for closest plane
    len = findIntersectionWithPlane(ray, vec3(1, 0, 0), pmax.x, out_intersect);
    if (len < minDist) { minDist = len; pClosest = out_intersect; }

    len = findIntersectionWithPlane(ray, vec3(-1, 0, 0), pmin.x, out_intersect);
    if (len < minDist) { minDist = len; pClosest = out_intersect; }

    len = findIntersectionWithPlane(ray, vec3(0, 1, 0), pmax.y, out_intersect);
    if (len < minDist) { minDist = len; pClosest = out_intersect; }

    len = findIntersectionWithPlane(ray, vec3(0, -1, 0), pmin.y, out_intersect);
    if (len < minDist) { minDist = len; pClosest = out_intersect; }

    len = findIntersectionWithPlane(ray, vec3(0, 0, 1), pmax.z, out_intersect);
    if (len < minDist) { minDist = len; pClosest = out_intersect; }

    len = findIntersectionWithPlane(ray, vec3(0, 0, -1), pmin.z, out_intersect);
    if (len < minDist) { minDist = len; pClosest = out_intersect; }

    out_intersect = pClosest;
    p = out_intersect.position;

    if (p.x >= pmin.x - EPS && p.x <= pmax.x + EPS && p.y >= pmin.y - EPS && 
        p.y <= pmax.y + EPS && p.z >= pmin.z - EPS && p.z <= pmax.z + EPS) return minDist;

    else return INFINITY;
}  

// Cylinder
float getIntersectOpenCylinder( Ray ray, vec3 center, vec3 axis, float len, float rad, out Intersection intersect ) {

    // calculate intersection point
    vec3 delP = ray.origin - center;

    vec3 temp1 = ray.direction - dot(ray.direction, axis) * axis;
    vec3 temp2 = delP - dot(delP, axis) * axis;

    float a = length(temp1) * length(temp1);
    float b = 2.0 * dot(temp1, temp2);
    float c = length(temp2) * length(temp2) - rad*rad;

    float t1 = (-b + sqrt(b*b-4.0*a*c))/(2.0*a);
    float t2 = (-b - sqrt(b*b-4.0*a*c))/(2.0*a);

    float t;
    if (t2 > EPS) t = t2;
    else if (t1 > EPS) t = t1;
    else return INFINITY;

    // cut to length
    intersect.position = rayGetOffset(ray, t);
    if (length(intersect.position - center)*length(intersect.position - center) - rad*rad > len*len)
        return INFINITY;

    // calculate intersection position and normal
    vec3 pDiff = intersect.position - center;
    float lenDiff = sqrt(length(pDiff)*length(pDiff)-rad*rad);
    intersect.normal = normalize(pDiff-normalize(axis)*lenDiff);
    return length(t * ray.direction);
}

float getIntersectDisc( Ray ray, vec3 center, vec3 norm, float rad, out Intersection intersect ) {
    float len = findIntersectionWithPlane(ray, norm, dot(norm, center), intersect);
    if (abs(length(intersect.position - center)) <= rad) return len;
    else return INFINITY;
}


float findIntersectionWithCylinder( Ray ray, vec3 center, vec3 apex, float radius, out Intersection out_intersect ) {
    vec3 axis = apex - center;
    float len = length( axis );
    axis = normalize( axis );

    Intersection intersect;
    float best_dist = INFINITY;
    float dist;

    // -- infinite cylinder
    dist = getIntersectOpenCylinder( ray, center, axis, len, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    // -- two caps
    dist = getIntersectDisc( ray, center, axis, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );
    dist = getIntersectDisc( ray,   apex, axis, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    return best_dist;
}
    
// Cone
float getIntersectOpenCone( Ray ray, vec3 apex, vec3 axis, float len, float radius, out Intersection intersect ) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 31 lines of code.
    vec3 center = axis*len + apex;

    vec3 v = normalize(ray.direction);
    vec3 va = normalize(axis);

    float alpha = atan(radius / len);
    vec3 delP = ray.origin - apex;

    vec3 temp1 = v - dot(v, va) * va;
    vec3 temp2 = delP - dot(delP, va) * va;
    float temp3 = dot(v, va);
    float temp4 = dot(delP, va);

    float a = cos(alpha)*cos(alpha) * length(temp1)*length(temp1) - sin(alpha)*sin(alpha)*temp3*temp3;
    float b = 2.0 * cos(alpha) * cos(alpha) * dot(temp1, temp2) - 2.0*sin(alpha)*sin(alpha)*temp3*temp4;
    float c = cos(alpha)*cos(alpha)*length(temp2)*length(temp2) - sin(alpha)*sin(alpha)*temp4*temp4;

    float t1 = (-b + sqrt(b*b-4.0*a*c))/(2.0*a);
    float t2 = (-b - sqrt(b*b-4.0*a*c))/(2.0*a);

    float t;
    if (t2 > 0.0) t = t2;
    else if (t1 > 0.0) t = t1;
    else return INFINITY;

    intersect.position = rayGetOffset(ray, t);
    vec3 pDiff = intersect.position - apex;

    if (length(pDiff)*length(pDiff) - radius*radius > len*len ||
         length(intersect.position - center) > len) return INFINITY;

    intersect.normal = normalize(pDiff - length(pDiff)/cos(alpha) * normalize(axis));
    return length(t * ray.direction);
}

float findIntersectionWithCone( Ray ray, vec3 center, vec3 apex, float radius, out Intersection out_intersect ) {
    vec3 axis   = center - apex;
    float len   = length( axis );
    axis = normalize( axis );
        
    // -- infinite cone
    Intersection intersect;
    float best_dist = INFINITY;
    float dist;

    // -- infinite cone
    dist = getIntersectOpenCone( ray, apex, axis, len, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    // -- caps
    dist = getIntersectDisc( ray, center, axis, radius, intersect );
    chooseCloserIntersection( dist, best_dist, intersect, out_intersect );

    return best_dist;
}

#define MAX_RECURSION 8

// simple random number generator http://www.ozone3d.net/blogs/lab/20110427/glsl-random-generator/
float rand(vec2 n) {
  return 0.5 + 0.5 * fract(sin(dot(n.xy, vec2(12.9898, 78.233)))* 43758.5453);
}

vec3 calculateSpecialDiffuseColor( Material mat, vec3 posIntersection, vec3 normalVector ) {
    if ( mat.special == CHECKERBOARD ) {
        
        // display square size similar to example
        posIntersection = posIntersection/8.0;

        if (mod(floor(posIntersection.x) + floor(posIntersection.y) + 
            floor(posIntersection.z), 2.0) == 0.0) return mat.color * 0.5;
        else return mat.color;
    }
    else if ( mat.special == MYSPECIAL ) {
        posIntersection = normalize(posIntersection);

        return mat.color;
        // return mat.color * (posIntersection.x + posIntersection.y + posIntersection.z);
    }

    return mat.color;
}

vec3 calculateDiffuseColor( Material mat, vec3 posIntersection, vec3 normalVector ) {
    // Special colors
    if ( mat.special != NONE ) {
        return calculateSpecialDiffuseColor( mat, posIntersection, normalVector ); 
    }
    return vec3( mat.color );
}

// check if position pos in in shadow with respect to a particular light.
// lightVec is the vector from that position to that light
bool pointInShadow( vec3 pos, vec3 lightVec ) {
    Ray ray;
    ray.origin = pos;
    ray.direction = normalize(lightVec);
    Material out_mat;
    Intersection out_intersect;

    float dist = rayIntersectScene( ray, out_mat, out_intersect );
    if (dist > EPS && dist < length(lightVec) + EPS) return true;
    return false;
}

// soft shadow
float pointShadowRatio( vec3 pos, vec3 lightVec ) {

    float count = 0.0;
    int k = 10;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            // randomly sample new light ray around original light
            float u = rand(vec2(-1.0, 1.0));
            float theta = rand(vec2(0.0, 2.0*PI));
            float x = sqrt(1.0 - u*u)*cos(theta);
            float y = sqrt(1.0 - u*u)*sin(theta);
            float z = u;
            vec3 newLighVec = vec3(x, y, z);
            if (!pointInShadow(pos, newLighVec)) {
                count += 1.0;
            }
        }
    }
    return count / float(k*k);
}

vec3 getLightContribution( Light light, Material mat, vec3 posIntersection, vec3 normalVector, vec3 eyeVector, bool phongOnly, vec3 diffuseColor ) {

    vec3 lightVector = light.position - posIntersection;
    
    if ( pointInShadow( posIntersection, lightVector ) ) {
        return vec3( 0.0, 0.0, 0.0 );
    }

    if ( mat.materialType == PHONGMATERIAL || mat.materialType == LAMBERTMATERIAL ) {
        vec3 contribution = vec3( 0.0, 0.0, 0.0 );

        // get light attenuation
        float dist = length( lightVector );
        float attenuation = light.attenuate * dist * dist;

        float diffuseIntensity = max( 0.0, dot( normalVector, lightVector ) ) * light.intensity;
        
        // glass and mirror objects have specular highlights but no diffuse lighting
        if ( !phongOnly ) {
            contribution += diffuseColor * diffuseIntensity * light.color / attenuation;
        }
        
        if ( mat.materialType == PHONGMATERIAL ) {
            vec3 refVector = reflect(lightVector, normalVector);
            float phongIntensity = pow(max( 0.0, dot( eyeVector, refVector ) ) * light.intensity, 0.1);
            // vec3 phongTerm = vec3( phongIntensity, phongIntensity, phongIntensity );
            vec3 phongTerm = vec3( 0.0, 0.0, 0.0 );

            contribution += phongTerm;
        }

        return contribution;
        // for soft shadow:
        return contribution * pointShadowRatio(posIntersection, lightVector);
    }
    else {
        return diffuseColor;
    }

}

vec3 calculateColor( Material mat, vec3 posIntersection, vec3 normalVector, vec3 eyeVector, bool phongOnly ) {
	vec3 diffuseColor = calculateDiffuseColor( mat, posIntersection, normalVector );

	vec3 outputColor = vec3( 0.0, 0.0, 0.0 ); // color defaults to black when there are no lights
	
    for ( int i=0; i<MAX_LIGHTS; i++ ) {

        if( i>=numLights ) break; // because GLSL will not allow looping to numLights
		
        outputColor += getLightContribution( lights[i], mat, posIntersection, normalVector, eyeVector, phongOnly, diffuseColor );
	}
	
	return outputColor;
}

// find reflection or refraction direction ( depending on material type )
vec3 calcReflectionVector( Material material, vec3 direction, vec3 normalVector, bool isInsideObj ) {
    if( material.materialReflectType == MIRRORREFLECT ) {
        return reflect( direction, normalVector );
    }
    // the material is not mirror, so it's glass.
    // compute the refraction direction...
    
    // ----------- STUDENT CODE BEGIN ------------
    // see lecture 13 slide ( lighting ) on Snell's law
    // the eta below is eta_i/eta_r
    float eta = ( isInsideObj ) ? 1.0/material.refractionRatio : material.refractionRatio;
    return refract( normalize(direction), normalize(normalVector), eta);

    // ----------- Our reference solution uses 11 lines of code.
    
    return reflect( direction, normalVector ); // return mirror direction so you can see something
    // ----------- STUDENT CODE END ------------
}

vec3 traceRay( Ray ray ) {
    Material hitMaterial;
    Intersection intersect;

    vec3 resColor  = vec3( 0.0, 0.0, 0.0 );
    vec3 resWeight = vec3( 1.0, 1.0, 1.0 );
    
    bool isInsideObj = false;

    for ( int depth = 0; depth < MAX_RECURSION; depth++ ) {
        
        float hit_length = rayIntersectScene( ray, hitMaterial, intersect );
            
        if ( hit_length < EPS || hit_length >= INFINITY ) break;

        vec3 posIntersection = intersect.position;
        vec3 normalVector    = intersect.normal;

        vec3 eyeVector = normalize( ray.origin - posIntersection );           
        if ( dot( eyeVector, normalVector ) < 0.0 )
            { normalVector = -normalVector; isInsideObj = true; }
        else isInsideObj = false;

        bool reflective = ( hitMaterial.materialReflectType == MIRRORREFLECT || 
                            hitMaterial.materialReflectType == GLASSREFLECT );
		vec3 outputColor = calculateColor( hitMaterial, posIntersection, normalVector, eyeVector, reflective );

        float reflectivity = hitMaterial.reflectivity;

        // check to see if material is reflective ( or refractive )
        if ( !reflective || reflectivity < EPS ) {
            resColor += resWeight * outputColor;
            break;
        }
        
        // bounce the ray
        vec3 reflectionVector = calcReflectionVector( hitMaterial, ray.direction, normalVector, isInsideObj );
        ray.origin = posIntersection;
        ray.direction = normalize( reflectionVector );

        // add in the color of the bounced ray
        resColor += resWeight * outputColor;
        resWeight *= reflectivity;
    }

    return resColor;
}

void main( ) {
    float cameraFOV = 0.8;
    vec3 direction = vec3( v_position.x * cameraFOV * width/height, v_position.y * cameraFOV, 1.0 );

    Ray ray;
	ray.origin    = vec3( uMVMatrix * vec4( camera, 1.0 ) );
    ray.direction = normalize( vec3( uMVMatrix * vec4( direction, 0.0 ) ) );

    // trace the ray for this pixel
    vec3 res = traceRay( ray );
    
    // paint the resulting color into this pixel
    gl_FragColor = vec4( res.x, res.y, res.z, 1.0 );
}

