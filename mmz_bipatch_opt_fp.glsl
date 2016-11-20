//------------------------------------------------------------------------
#include "data/shaders/heightmap/base.glsl"

// -------------------------------------------------------
uniform float lodBaseDistance;

// -------------------------------------------------------
int lodLevel = -1;

#if 0
//-----------------------------------------------------------------------------------------
real _computeUFromV(real v, real2 a12, real2 b12, real2 c12, real2 d12)
{
    return (v * (c12.x - c12.y) + (d12.x - d12.y)) / (v * (a12.y - a12.x) + (b12.y - b12.x));

    /*real a = v * a12.y + b12.y;
    real b = v * (a12.y - a12.x) + b12.y - b12.x;
    if (abs(b) >= abs(a))
        return (v * (c12.x - c12.y) + d12.x - d12.y) / b;
    else
        return (-v * c12.y - d12.y) / a;*/
}

//-----------------------------------------------------------------------------------------
real2 _compute2UFrom2V(real2 vv, real2 a12, real2 b12, real2 c12, real2 d12)
{
    real2 aa = vv * a12.y + b12.y;
    real2 bb = vv * (a12.y - a12.x) + b12.y - b12.x;
    //return (vv * (c12.x - c12.y) + d12.x - d12.y) / bb;
    return -(vv * c12.y + d12.y) / (vv * a12.y + b12.y);
}

//-----------------------------------------------------------------------------------------
real3 _computePuv(real3 r, real3 q, real2 uv, real3 a, real3 b, real3 c, real3 d)
{
    return uv.s * uv.t * a + uv.s * b + uv.t * c + d;
}

//-----------------------------------------------------------------------------------------
real2 _compute2Intersections(real3 r, real3 q, real4 u, real3 a, real3 b, real3 c, real3 d)
{
    if (abs(q.x) > 0)
        return (u.xz * u.yw * a.xx + u.xz * b.xx + u.yw * c.xx + d.xx - r.xx) / q.x;
    if (abs(q.y) > 0)
        return (u.xz * u.yw * a.yy + u.xz * b.yy + u.yw * c.yy + d.yy - r.yy) / q.y;
    if (abs(q.z) > 0)
        return (u.xz * u.yw * a.zz + u.xz * b.zz + u.yw * c.zz + d.zz - r.zz) / q.z;
}


//-----------------------------------------------------------------------------------------
real4 rayBilinearPlaneIntersectionOptimized(real3 r, real3 q, real3 p00, real3 p01, real3 p10, real3 p11)
{
    // precompute some values for substitution
    real3 a = p11 - p10 - p01 + p00;
    real3 b = p10 - p00;
    real3 c = p01 - p00;
    real3 d = p00;

    // precompute other variables for substitution
    real2 a12 = a.xy * q.z - a.z * q.xy;
    real2 b12 = b.xy * q.z - b.z * q.xy;
    real2 c12 = c.xy * q.z - c.z * q.xy;
    real2 d12 = (d.xy - r.xy) * q.z - (d.z - r.z) * q.xy;
    
    // precompute a,b,c terms of an quadratic equation
    real _a = a12.y * c12.x - a12.x * c12.y;
    real _b = a12.y * d12.x - a12.x * d12.y + b12.y * c12.x - b12.x * c12.y;
    real _c = b12.y * d12.x - b12.x * d12.y;
    
    // special case if _a == 0, hence the bilinear patch is planar
    if (abs(_a) == 0)// < 0.0000001)
    {
        // compute uv
        real2 uv;
        uv.t = -_c / _b;
        uv.s = _computeUFromV(uv.t, a12, b12, c12, d12);

        // check if values are valid and there is an intersection
        if (any(lessThan(uv, 0)) || any(greaterThan(uv, 1))) return real4(r, 0);
        
        // compute intersection point
        return real4(_computePuv(r, q, uv, a, b, c, d), 1);
    }

    // compute discreminant of the quadratic equation
    real _d = _b * _b - 4.0 * _a * _c;
    //real _d = _b * _b / (4.0 * _a * _a) - _c / _a;

    // if there is no intersection at all
    if (_d < 0)
    {
        return real4(r, 0);

    // if there is only one result, then 
    /*}else if (_d == 0)
    {
        // uv
        real2 uv;

        // compute uv
        uv.t = -0.5 * _b / _a;
        uv.s = _computeUFromV(uv.t, a12, b12, c12, d12);

        // check if values are valid and there is an intersection
        if (any(lessThan(uv, 0)) || any(greaterThan(uv, 1))) return real4(r, 0);

        // compute intersection
        return real4(_computePuv(r, q, uv, a, b, c, d), 1);

    // if there exists both intersections*/
    }else
    {
        // precompute some values
        real _t = 0.5 * sqrt(_d) / _a;
        real _s = -0.5 * _b / _a;
        int lp0 = 1;
        int lp1 = 1;
    
        // compute uv for both intersections at the same time 
        real4 uvuv;
        uvuv.yw = real2(_s - _t, _s + _t);
        uvuv.xz = _compute2UFrom2V(uvuv.yw, a12, b12, c12, d12);
    
        // if uv coordinates are not valid, then do not compute intersection point
        if (any(lessThan(uvuv.xy, 0)) || any(greaterThan(uvuv.xy, 1))) lp0 = 0;
        if (any(lessThan(uvuv.zw, 0)) || any(greaterThan(uvuv.zw, 1))) lp1 = 0;

        if (lp0 == 0 && lp1 == 0) return real4(r, 0);
    
        // compute both intersection of the ray and the patch
        real2 tt = _compute2Intersections(r, q, uvuv, a, b, c, d);
        real3 p0 = r + tt.x * q;
        real3 p1 = r + tt.y * q;

        // get the nearest allowed intersection
        if (lp0 == 1 && lp1 == 1)
            if (tt.x < tt.y) return real4(p0, 1);
            else return real4(p1, 1);
        else if (lp0 == 0 && lp1 == 1) return real4(p1,1);
        else if (lp0 == 1 && lp1 == 0) return real4(p0,1);
            
    
        // compute the nearest intersection in the positive direction of the ray
        /*real _s = -0.5 * _b / _a;

        // compute uv coordinates
        real2 uv;
        uv.t = _s - sign(_s) * sqrt(_d);
        uv.s = _computeUFromV(uv.t, a12, b12, c12, d12);
        real3 p0 = _computePuv(r, q, uv, a, b, c, d);

        // if uv coordinates are not valid, then do not compute intersection point
        if (any(lessThan(uv.st, 0)) || any(greaterThan(uv.st, 1))) return real4(r, 0);
                
        return real4(p0, 1);*/
    }
}

/**
 * Same as computeRayPatchIntersection, but use the fact, that our
 * texture is linearly interpolated
 **/
real4 computeRayPatchIntersectionLinear(real3 r0, real3 rD, real4 _value, vec4 Pd, float nextIntersection, real2 sizeBig)
{
    // setup the four corners of the bilinear patch 
    real3 p[4];

    p[0].xz = Pd.xy;
    p[1].xz = Pd.xw;
    p[2].xz = Pd.zy;
    p[3].xz = Pd.zw;

    p[0].y = _value.x;
    p[1].y = _value.z;
    p[2].y = _value.y;
    p[3].y = _value.w;

    // compute rays intersection with the bilinear patch
    return rayBilinearPlaneIntersectionOptimized(r0.xyz, rD, p[0], p[1], p[2], p[3]);
}

#else

real4 computeRayPatchIntersectionLinear(real3 r0, real3 rD, real4 _value, vec4 Pd, float nextIntersection, real2 sizeBig)
{
    // here we store the result
    real4 result = 0;

    // linear search steps
    const int LINEAR_STEPS = 5;
    float ds = nextIntersection / LINEAR_STEPS;

    // perform linear search of the intersection point
    for (int i=0; i < LINEAR_STEPS; i++)
    {
        // move ray 
        r0 = r0 + ds * rD;

        // get height on the current position
        //real4 Texel = texture2DLod(heightmapTex, r0.xz, 0).xyzw * maxHeight;
        //real Zmin = max(max(Texel.x, Texel.y), max(Texel.z, Texel.w));
        real Zmin = texture2DLod(heightmapTex, r0.xz + 0.0 / sizeBig.xy, 0.0).x * maxHeight;
        //real Zmin = max(max(Texel.x, Texel.y), max(Texel.z, Texel.w));

        // check if we found intersection point 
        if (Zmin >= r0.y)
        {
            result.w = 1;
            result.xyz = r0;
            break;
        }         
        iterationCounter ++;
    }

    // size of the search window
    float size = ds;
    float step = 1.0;

    // now do a binary search to find the correct intersection point
    for (int i=0; i < 3; i++)
    {
        size *= 0.5;

        // get current texel value 
        //real4 Texel = texture2DLod(heightmapTex, r0.xz, 0).xyzw * maxHeight;
        //real Zmin = max(max(Texel.x, Texel.y), max(Texel.z, Texel.w));
        real Zmin = texture2DLod(heightmapTex, r0.xz + 0.0 / sizeBig.xy, 0.0).x * maxHeight;

        // check if we allowed to move to the boundary
        if (r0.y < Zmin){
            step = -1.0;
        }else{
            step = 1.0;
        }

        // move the ray on the current depth
        r0 = r0 + step * size * rD;

        iterationCounter ++;
    }
    result.xyz = r0;
    return result;
}

#endif

/**
 * Compute ray propagation until it stops and returns ray position
 **/
vec4 propagateRay(vec3 _r0, vec3 rS, vec3 rD, vec3 eyePosLocal, vec3 _min = vec3(0,0,0), vec3 _max = vec3(1,maxHeight,1))
{
    // we start from this level 
    int stepCounter = 0;
    real4 r0 = real4(_r0, 0);

    // compute starting level
    real2 texSizeBig = textureSize2D(heightmapTex, 0);
    int level = log2(texSizeBig.x);

    // texture and inverse texture size
    real2 texSize = real2(1,1);

    // here we store ray coordinates in texture size domain
    real2 rC = r0.xz * texSize;
    real srDy = sign(rD.y);

    // repeat
    while(stepCounter++ < MAX_STEP)
    {
        // check if pos is on horizontal and/or vertical plane
        //vec2 hv = frac(rC) < 0.000001;
        real2 hv = floor(1 - frac(rC));
    
        // compue sampling coordinates based on this information
        ivec2 st = (floor(rC)  + floor(min(rD.xz, vec2(0,0))) * hv);
 
        // sample the sampler and return the value 
        real4 Texel = texelFetch2D(heightmapTex, st, max(level, 0)).xyzw * maxHeight;

        // set size of the heightmap below the ray
        real Zmin = max(max(Texel.x, Texel.y), max(Texel.z, Texel.w));

        // compute distancies of mipmap boundaries
        real4 Pd;
        Pd.xy = real2(st.st + 0) / real2(texSize.xy);
        Pd.zw = real2(st.st + 1) / real2(texSize.xy);

        // we store here the nearest intersection
        real intersection;
        real nextIntersection;

        // decide on the direction of the ray which intersection to compute
        if (rD.x < 0)
        {
            if (rD.z < 0)
            {
                // there are two possible intersection (left and top)
                real2 t = (Pd.xy - r0.xz) / rD.xz;
                intersection = t.x;
                intersection = min(intersection, t.y);
            }else
            {
                // there are two possible intersection (left and bottom)
                real2 t = (Pd.xw - r0.xz) / rD.xz;
                intersection = t.x;
                intersection = min(intersection, t.y);
            }
        }else
        {
            if (rD.z < 0)
            {
                // there are two possible intersection (right and top)
                real2 t = (Pd.zy - r0.xz) / rD.xz;
                intersection = t.x;
                intersection = min(intersection, t.y);
            }else
            {
                // there are two possible intersection (right and bottom)
                real2 t = (Pd.zw - r0.xz) / rD.xz;
                intersection = t.x;
                intersection = min(intersection, t.y);
            }
        }

        // compute intersection with the Z-plane
        nextIntersection = intersection;
        intersection = min(srDy * (r0.y - Zmin) / rD.y, intersection);

        
        // based on the distance of the ray to the eye, we compute the lod level 
        //if (dot(r0.xyz - eyePosLocal, r0.xyz - eyePosLocal) > pow(lodBaseDistance, 2))
        if (distance(r0.xyz, eyePosLocal) > lodBaseDistance)
        {
            lodLevel ++;
            lodBaseDistance *= 2.0; 
        }

        // stop iteration only if ray falls below of the z map
        if (level == lodLevel && r0.y <= Zmin )
        {
            // do only compute the intersection, if lod level is 0
            if (lodLevel >= 0){r0.w = 1; break; }
            else 
            {
                // now compute the interpolated bilinar patch corners
                real4 pos = computeRayPatchIntersectionLinear(r0, rD, Texel, Pd, nextIntersection, texSizeBig);
    
                // if intersection was successfull, then move the ray there and stop 
                if (pos.w > 0)
                {
                    r0 = pos;
                    break;
                }else{
                    // move the ray on the next best position and continue the propagation
                    intersection = nextIntersection;
                }
            }
        }

        // we only move the ray if it intersects with the mipmaps boundaries or in level 0
        if ((r0.y > Zmin ) || level == lodLevel)
        {           
            // move ray to the nearest intersection point
            r0.xyz = r0.xyz + intersection * rD;

            // compute new coordinates
            real2 _r = r0.xz * real2(texSize);
        
            // check if divided by two we get a rest
            _r = _r * 0.5;
            _r = floor(1 - frac(_r));
            real inc = max(_r.x, _r.y);
            level += inc;

        }else{

            // go deeper in the mipmap
            level = max(level - 1, lodLevel);            
        }

        // else decrement level
        level = max(level, lodLevel);

        // get current texture size
        texSize = textureSize2D(heightmapTex, max(level, 0));
        rC = r0.xz * texSize;

        // iteration counter
        iterationCounter ++;

        // check if we are moving outside
        if (any(greaterThanEqual(vec2(rC), vec2(texSize)))) break;
        if (any(lessThanEqual(vec3(r0),0.0))) break;    
    }

    return r0;
}
