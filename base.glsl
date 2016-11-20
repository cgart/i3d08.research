#ifndef HM_BASE_GLSLS
#define HM_BASE_GLSLS

// -------------------------------------------------------
#include "data/shaders/utils.glsl"

// -------------------------------------------------------
#define HeightMapSampler sampler2D

// -------------------------------------------------------
uniform HeightMapSampler heightmapTex;
uniform sampler2D heightmapDecalTex;
uniform sampler2D heightmapNormalTex;

// -------------------------------------------------------
typedef half real;
typedef vec2 real2;
typedef vec3 real3;
typedef vec4 real4;

// -------------------------------------------------------
#define LEFT 0
#define TOP 1
#define RIGHT 2
#define BOTTOM 3
#define ZPLANE 4
#define PLANE_COUNT 5
#define MAX_STEP 4128

// -------------------------------------------------------
int iterationCounter = 0;

/**
 * Compute largest mipmap level which can be accessed from the current position
 **/
inline int computeLargestMipmapLevel(int fish)
{
    return 1 - fish & 1;
    //if (fish == 0) return 0;
    /*int counter = 0;
    while(fish && (fish & 1 == 0))
    {
        //result = result + (fish & 1 == 0) ? 1 : 0;
        fish >>= 1;
        counter ++;
    }
    return counter;*/
}

//! Compute intersection between ray and a plane
inline float rayPlaneIntersection(vec3 r0, vec3 rD, vec3 Pn, float Pd)
{    
    // check if parallel
    //float q = dot(rD, Pn);
    
    // compute intersection
    //float t = -(dot(r0, Pn) + Pd) / dot(rD, Pn);

    return -(dot(r0, Pn) + Pd) / dot(rD, Pn);
}


//! Sample a sampler based on the direction
vec4 texelFetch2DByDirection(HeightMapSampler tex, vec2 pos, vec2 dir, vec2 texSize, int level, out ivec2 st)
{
    // convert pos into texture size domain
    vec2 sst = pos * texSize;

    // check if pos is on horizontal and/or vertical plane
    //vec2 hv = frac(sst) < 0.000001;
    vec2 hv = floor(1 - frac(sst));

    // compue sampling coordinates based on this information
    st = (sst + floor(min(dir.xy, vec2(0,0))) * hv);//  * invTexSize;

    // sample the sampler and return the value 
    //return texelFetch3D(tex, ivec3(st,level), 0 ).x;//yzw;
    return texelFetch2D(tex, st, level).xyzw;
    //return texture2DLod(tex, (st + 0.5) / texSize, level).xyzw;
}

//! Compute intersection and thus starting and endposition of the ray and bounding box of the volume
inline vec4 computeRayHeightmapAABBIntersection(vec3 eyePos, vec3 rD, out vec3 rS, vec3 _min = 0, vec3 _max = vec3(1,maxHeight,1))
{
    // our starting position
    vec4 r0 = 0;

    // compute intersection with the bounding box of the heightfield
    vec3 iner = rayBoxIntersection(eyePos, rD, _min, _max);
    r0.w = iner.z;

    // if our starting points is not in the box, then compute first and last position
    if (any(lessThan(eyePos, _min)) || any(greaterThan(eyePos, _max)))
    {
        r0.xyz = eyePos + (iner.x + 0.0001) * rD;
        rS = eyePos + (iner.y - 0.0001) * rD;        
    }else{
        r0.xyz = eyePos;
        rS = eyePos + (iner.x - 0.0001) * rD;
    }

    return r0;
}

#endif
