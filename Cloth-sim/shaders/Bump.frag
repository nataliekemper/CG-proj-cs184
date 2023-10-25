#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_2;
uniform vec2 u_texture_2_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
  // You may want to use this helper function...
  return texture(u_texture_2, uv).r;
}

void main() {
  // Compute displaced normal
  vec4 b = vec4(cross(v_tangent.xyz, v_normal.xyz), 0);

  float dU = (h(vec2(v_uv.x + 1 / u_texture_2_size.x, v_uv.y)) - h(vec2(v_uv.x, v_uv.y))) * u_normal_scaling * u_height_scaling;
  float dV = (h(vec2(v_uv.x, v_uv.y + 1 / u_texture_2_size.y)) - h(vec2(v_uv.x, v_uv.y))) * u_normal_scaling * u_height_scaling;

  vec3 n0 = vec3(-dU, -dV, 1);
  vec4 n = v_tangent * n0.x + b * n0.y + v_normal * n0.z;

  // (Placeholder code. You will want to replace it.)
  // parameters
  vec4 kd = u_color;
  float ka = 0.1;
  float ks = 0.5;
  vec4 Ia = vec4(1, 1, 1, 1);
  float p = 100;
  
  vec4 r = vec4(u_light_pos, 0.0) - v_position;
  vec4 I = vec4(u_light_intensity, 0.0);
  vec4 l = r / length(r.xyz);
  vec4 v = vec4(u_cam_pos, 0.0) - v_position;
  v = v / length(v.xyz);
  vec4 h = (v + l) / length((v + l).xyz);

  vec4 ambient = ka * Ia;
  vec4 diffuse = kd * I / dot(r, r) * max(0.0, dot(n, l));
  vec4 specular = ks * I / dot(r, r) * pow(max(0.0, dot(n, h)), p);
  out_color = ambient + diffuse + specular;
  out_color.a = 1;
}
