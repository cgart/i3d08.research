
#ifndef UTILS_GLSL
#define UTILS_GLSL

/**
* Just as reminder. you should rather use uniform values like
* inFragDC1=(far * near / (near - far));
* inFragDC2=far / (far - near);
* and go for a simple gl_FragDepth = inFragDC2 + inFragDC1/eyeCoordZ;
*/
float eyeZToFragDepth(float z, float near, float far){
    return ((-far / (far - near) * z - far * near / (near - far)) / -z);
    //return far / (far - near) + ((far * near / (near - far)) / z);
    //return (far + near) / (far - near) + ((-2.0 * far * near / (far - near)) / z);
}

float selectMipmapLevel(vec2 uv, vec2 textureSize)
{
    vec2 dx = ddx(uv * textureSize.x);
    vec2 dy = ddy(uv * textureSize.y);
    float d = max( dot(dx, dx), dot(dy, dy) );

    return log2( sqrt(d) );
}

/**
 * Compute fresnel term based on the direction and the normal 
 * Assume given parameters are normalized.
 **/
vec2 computeFresnelTerm(vec3 normal, vec3 direction)
{
    vec2 RT;
    
    RT.x = 1 / pow(1 + dot(normal, direction), 2.0);
    RT.y = (1 - RT.x);

    return RT;
}

//! Convert normal map coordinates
vec3 textureNormalToWorldNormal (vec3 v)
{
    return (v.xzy - 0.5)*2.0;
}

//! Compute Diffuse Shading
vec3 computeDiffuse(in vec3 N, in vec3 light, in vec3 v, in vec3 _min)
{
    vec3 diffuseL = max(dot(N, normalize(light - v) ), _min);

    return diffuseL;
}

/**
 * Ray-box intersection using IEEE numerical properties to ensure that the
 * test is both robust and efficient, as described in:
 *
 *      Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
 *      "An Efficient and Robust Ray-Box Intersection Algorithm"
 *      Journal of graphics tools, 10(1):49-54, 2005
 *
 * @return vec3(tmin, tmax, 0 or 1 if there was a cut or not)
 *
 */
vec3 rayBoxIntersection(vec3 r0, vec3 rD, vec3 _min, vec3 _max)
{
    float tmin, tmax, tymin, tymax, tzmin, tzmax;

    // inverse direction to catch float problems
    vec3 invrd = 1.0 / rD;

    if (invrd.x >= 0)
    {
        tmin = (_min.x - r0.x) * invrd.x;
        tmax = (_max.x - r0.x) * invrd.x;
    }else
    {
        tmin = (_max.x - r0.x) * invrd.x;
        tmax = (_min.x - r0.x) * invrd.x;
    }

    if (invrd.y >= 0)
    {
        tymin = (_min.y - r0.y) * invrd.y;
        tymax = (_max.y - r0.y) * invrd.y;
    }else
    {
        tymin = (_max.y - r0.y) * invrd.y;
        tymax = (_min.y - r0.y) * invrd.y;
    }

    if ( (tmin > tymax) || (tymin > tmax) ) return 0;
    if ( tymin > tmin) tmin = tymin;
    if ( tymax < tmax) tmax = tymax;

    
    if (invrd.z >= 0)
    {
        tzmin = (_min.z - r0.z) * invrd.z;
        tzmax = (_max.z - r0.z) * invrd.z;
    }else
    {
        tzmin = (_max.z - r0.z) * invrd.z;
        tzmax = (_min.z - r0.z) * invrd.z;
    }

    if ( (tmin > tzmax) || (tzmin > tmax) ) return 0;
    if ( tzmin > tmin) tmin = tzmin;
    if ( tzmax < tmax) tmax = tzmax;

    // check if values are valid
    if (tmin < 0) tmin = tmax;
    if (tmax < 0) return 0;

    return vec3(tmin, tmax, 1);
}

#endif

