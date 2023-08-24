#version 450
#extension GL_EXT_samplerless_texture_functions : enable
layout(set = 0, binding = 0) uniform sampler samp;
layout(set = 0, binding = 1) uniform texture2D src;

layout(location = 0) in vec2 fragUv;

layout(location = 0) out vec4 outColor;
float gaussianKernel[5][5] = {
    { 1,  4,  6,  4, 1 },
    { 4, 16, 24, 16, 4 },
    { 6, 24, 36, 24, 6 },
    { 4, 16, 24, 16, 4 },
    { 1,  4,  6,  4, 1 }
};

vec3 bilateralFilter(vec2 uv) {
    vec3 colorCenter=texture(sampler2D(src,samp),uv).xyz;

    vec2 textureSize = textureSize(src, 0);
    float offsetX=1.0/textureSize.x;
    float offsetY=1.0/textureSize.y;
    
    vec3 colorSum=vec3(0);
    float weightSum=0;
    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++){
            vec2 uvOff=uv+vec2((i-2)*offsetX,(j-2)*offsetY);
            vec3 color=texture(sampler2D(src,samp),uvOff).xyz;
            
            vec3 test=colorCenter-color;
            test*=test;
            if(test.x>0.2&&test.y>0.2&&test.z>0.2){
                continue;
            }
            
            colorSum+=texture(sampler2D(src,samp),uvOff).xyz*gaussianKernel[i][j];
            weightSum+=gaussianKernel[i][j];

        }
    }
    return colorSum/weightSum;
    
}

void main() {

    outColor = vec4(bilateralFilter(fragUv),1.0);
}