
//--------------------------------------------------
varying vec3 Normal;
varying vec4 Vertex;
varying vec4 Color;
varying vec4 HPos;
varying vec3 Tangent;
varying vec3 Binormal;
varying vec3 eyeVec;
varying vec3 tanEyeVec;
varying vec3 tanEyePos;
varying vec3 tanLightVec;
varying vec3 tanHalfViewVec;

//--------------------------------------------------
uniform mat4 localToWorld;
uniform mat4 worldToLocal;
uniform mat4 localToWorldIT;
uniform mat4 worldToLocalIT;
uniform vec3 eyePos;
uniform float maxHeight;
uniform vec3 sunPosition;

// -------------------------------------------------------
attribute vec3 vertexTangent;



/**
 * Compute the shading and bypass some data to the fragment program
 **/
void main(void)
{
    // bypass the data
    gl_TexCoord[0].xy = gl_MultiTexCoord0.xy;
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    HPos = gl_Position;
    Color = gl_Color;

    // vertex in object space
    Vertex = gl_Vertex;

    // compute TBN matrix (object space to tangent space)
    Normal = normalize(gl_Normal);
    Tangent = normalize(vertexTangent);
    Binormal = normalize(cross(Normal, Tangent));
    mat3 TBN = mat3(Tangent, Normal, Binormal);
    TBN = transpose(TBN);

    // transform eye position to object space
    vec4 _eye = gl_ModelViewMatrixInverse * vec4(0,0,0,1);
    vec3 eye = _eye.xyz / _eye.w;
    vec3 eyeVec = normalize(Vertex.xyz - eye);

    // transform light vector into object space
    vec4 _sunPosition = vec4(sunPosition, 1);
    _sunPosition = worldToLocal * _sunPosition; // NOTE: maybe wrong for reflection
    sunPosition.xyz = _sunPosition.xyz / _sunPosition.w;
    vec3 lightDir = normalize(sunPosition - Vertex.xyz);

    // compute eye position in tangent space
    //tanEyePos = TBN * eye;

    // compute view direction in tangent space
    //tanEyeVec = TBN * eyeVec;

    // transform vertex to tangent space
    //Vertex.xyz = TBN * Vertex.xyz;


    // transform light vector into tangent space
    tanLightVec = TBN * lightDir;

    // compute half view vector (for specular lighting)
    tanHalfViewVec = TBN * normalize(lightDir - eyeVec);

    // transform eye pos into object space
    //vec4 _eye = gl_ModelViewMatrixInverse * vec4(0,0,0,1);
    //vec3 eye = _eye.xyz / _eye.w;
    //tanEyePos = eye;// + vec3(0.5, 0, 0.5);

    // compute view vector in object space
    //eyeVec = normalize(Vertex - eye);
    //tanEyePos = eye + vec3(0.5, 0, 0.5);

    //Vertex.xz = gl_TexCoord[0].xy;
    //Vertex.y = maxHeight;

    // transform view vector from object space into tangent space
    //tanEyeVec = eyeVec * TBN;
    //tanEyeVec = normalize(tanEyeVec);

    // compute eye vector in tangent space
    //tanEyeVec = normalize(TBN * eyeDirection);

    // compute eye position in tangent space
    //tanEyePos = Vertex - eyeDirection;
    //tanEyePos *= TBN;
}

