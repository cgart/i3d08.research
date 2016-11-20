// Heightmap texture used to create mipmapping
uniform sampler2D heightmapTex;

// Our output texture is here as input, but we read inputs from the levels below
uniform sampler2D outputTex;

// current mipmap level where we render the output
uniform float osgppu_MipmapLevel;

uniform float time;
uniform float heightCorrection;

uniform float osgppu_ViewportWidth;
uniform float osgppu_ViewportHeight;

//varying out uvec4 fragData;

/*
 * Compute minimum mipmap for a given 2D texture.
 */
void main(void)
{
    // if we are in the mipmaplevel 0, then jsut copy the values from the input texture
    if (osgppu_MipmapLevel < 0.5)
    {
        // get size of the input texture 
        ivec2 size = ivec2((int)osgppu_ViewportWidth, (int)osgppu_ViewportHeight);//textureSize2D(outputTex, 0);
        ivec2 oCoord = ivec2(gl_TexCoord[0] * size);

        // compute input coordinates on which we will sample 
        ivec2 iCoord[4];

        iCoord[0] = oCoord + ivec2(0,0);
        iCoord[1] = oCoord + ivec2(1,0);
        iCoord[2] = oCoord + ivec2(0,1);
        iCoord[3] = oCoord + ivec2(1,1);
        
        // sample the input texture on that positions
        float c[4];
        for (int i=0; i < 4; i++)
        {
            iCoord[i] = clamp(iCoord[i], 0, size - 1);
            c[i] = texelFetch2D(heightmapTex, iCoord[i], 0) + heightCorrection;
        }

        // store the results into the level 0 of our output texture
        gl_FragData[0].xyzw = vec4(c[0], c[1], c[2], c[3]);
       
        // compute maximum of all the values
        //float m = max( max(c[0], c[1]), max(c[2], c[3]) );
    
        // compute result for the according layer
        //gl_FragData[0].xyzw = m;//max( max(m.x, m.y), max(m.z, m.w) );

    // for other levels we do compute the maximal value of the corner points stored in channels
    }else{

        // get texture sizes of the previous level
        ivec2 size = textureSize2D(outputTex, (int)osgppu_MipmapLevel - 1);
        //ivec2 size = ivec2((int)osgppu_ViewportWidth, (int)osgppu_ViewportHeight);//textureSize2D(outputTex, 0);
        //size.st *= 2;

        // this our starting sampling coordinate 
        ivec2 iCoord = ivec2(gl_TexCoord[0] * size);
    
        // create offset for the texel sampling (TODO check why -1 seems to be correct)
        ivec2 st[4];
        st[0] = iCoord + ivec2(0,0);
        st[1] = iCoord + ivec2(-1,0);
        st[2] = iCoord + ivec2(0,-1);
        st[3] = iCoord + ivec2(-1,-1);
        
        // retrieve 4 texels from the previous mipmap level
        vec4 c[4];
        for (int i=0; i < 4; i++)
        {
            // map texels coordinates, such that they do stay in defined space
            st[i] = clamp(st[i], 0, size.xy - 1);
            
            // get texel from the previous mipmap level
            c[i] = texelFetch2D(outputTex, st[i], (int)osgppu_MipmapLevel - 1).xyzw;
        }
        
        // compute maximum of all the values
        vec4 m = max( max(c[0].xyzw, c[1].xyzw), max(c[2].xyzw, c[3].xyzw) );
    
        // compute result for the according layer
        gl_FragData[0].xyzw = max( max(m.x, m.y), max(m.z, m.w) );
    }
}
