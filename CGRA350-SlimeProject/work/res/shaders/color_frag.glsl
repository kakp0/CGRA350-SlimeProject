// res/shaders/color_frag.glsl

#version 330 core

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec4 uColor; // CHANGED: from vec3 to vec4 to include alpha

// viewspace data (this must match the output of the fragment shader)
in VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} f_in;

// framebuffer output
out vec4 fb_color;

void main() {
	// calculate lighting (hack)
	vec3 eye = normalize(-f_in.position);
	float light = abs(dot(normalize(f_in.normal), eye));
	
    // CHANGED: Lighting is now applied to the .rgb part of the color
	vec3 color = mix(uColor.rgb / 4.0, uColor.rgb, light);
	
	// output to the frambuffer
    // CHANGED: The alpha value now comes from uColor.a instead of being hardcoded to 1.0
	fb_color = vec4(color, uColor.a);
}