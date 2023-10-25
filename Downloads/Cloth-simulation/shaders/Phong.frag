#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;

out vec4 out_color;

void main() {

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
  vec4 n = v_normal;
  vec4 h = (v + l) / length((v + l).xyz);

  vec4 ambient = ka * Ia;
  vec4 diffuse = kd * I / dot(r, r) * max(0.0, dot(n, l));
  vec4 specular = ks * I / dot(r, r) * pow(max(0.0, dot(n, h)), p);
  out_color = ambient + diffuse + specular;

  out_color.a = 1;
}
