// -------------------------------------------------------
uniform sampler2D gradTex;

//--------------------------------------------------
uniform mat4 localToWorld;
uniform mat4 worldToLocal;
uniform mat4 localToWorldIT;
uniform mat4 worldToLocalIT;
uniform float maxHeight;
uniform float nearPlane;
uniform float farPlane;
uniform vec3 eyePos;

uniform float colorFactor;

// -------------------------------------------------------
varying vec3 Normal;
varying vec4 Vertex;
varying vec4 Color;
varying vec4 HPos;
varying vec3 Tangent;
varying vec3 Binormal;
varying vec3 tanEyePos;
varying vec3 tanEyeVec;

// -------------------------------------------------------
#include "data/shaders/utils.glsl"
#include "data/shaders/heightmap/mmz_bipatch_opt_fp.glsl"
//#include "data/shaders/heightmap/linear_binary.glsl"


// -------------------------------------------------------
void main(void)
{
    // compute ray direction
    vec3 r0 = Vertex;
    vec4 _eP = gl_ModelViewMatrixInverse * vec4(0,0,0,1);
    vec3 eP = _eP.xyz / _eP.w;

    vec3 rD = normalize(r0 - eP);
    vec3 rS = r0;
    eP.xz += 0.5;

    // compute starting and end position of the ray
    vec4 sr0 = computeRayHeightmapAABBIntersection(eP, rD, rS);

    // continue only on valid intersection
    if (sr0.w > 0)
    {
        // propagate the ray until it stops
        vec4 rr0 = propagateRay(sr0.xyz, rS, rD, eP);

        //float intersection = sr0.y * (rS.y - maxHeight) / rD.y;
        //vec4 rr0 = vec4(rS + rD * intersection, 1);

        // compute color of the resulting ray position
        Color.rgba = texture2D(heightmapDecalTex, rr0.xz, 0.0).rgba;
        Normal = texture2D(heightmapNormalTex, rr0.xz, 0.0).xzy * 2.0 - 1.0;
        Normal = normalize(Normal);

        // compute lighting for both lights
        vec3 diffuse = 0.0;
        for (int i=0; i < 2; i++)
        {
            // opengl lights are in view space, so convert them to local space
            vec4 lpos = gl_ModelViewMatrixInverse * gl_LightSource[i].position;
            lpos.xyz /= lpos.w;
            lpos.xz += 0.5;

            // compute diffuse lighting 
            vec3 lrD = lpos.xyz - rr0.xyz;
            lrD = normalize(lrD);
            diffuse += max(dot(Normal, lrD), vec3(0,0,0)) * gl_LightSource[i].diffuse.rgb;
        }

        // setup terrain color
        gl_FragColor.rgb = diffuse * Color.rgb * colorFactor ;
        gl_FragColor.a = rr0.w;

        // compute the depth of the pixel
        rr0.xz -= 0.5;
        vec4 _r0 = gl_ModelViewProjectionMatrix * vec4(rr0.xyz , 1);
        gl_FragDepth = eyeZToFragDepth (_r0.z , nearPlane, farPlane);
    }
}
