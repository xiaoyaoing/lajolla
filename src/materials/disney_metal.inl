#include "../microfacet.h"

Real get_D(Vector3 h,Vector2 alpha){
    Real alphaX = alpha.x;
    Real alphaY = alpha.y;
    Real ax2 = alphaX * alphaX;
    Real ay2 = alphaY * alphaY;
    Real D = M_PI * alphaX *alphaY * pow(h.x/ax2+h.y/ay2+h.z,2);
    D = 1/D;
    return D;
}

Vector2 get_alpha(Real roughness,Real anis){
    Real aspect = sqrt(1-0.9*anis);
    Real alphaX = max(0.0001,pow(roughness,2)/aspect);
    Real alphaY = max(0.0001,pow(roughness,2)*aspect);
    return {alphaX,alphaY};
}

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0 ||
            dot(vertex.geometry_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Vector3  wh = normalize((dir_in + dir_out));
    Spectrum baseColor = eval(bsdf.base_color,vertex.uv,0,texture_pool);

    //Fresnel item
    Spectrum F = baseColor + (1-baseColor) * pow((1- absDot(wh,frame.n)),5);
    //NDF item
    Real anisotropic = eval(bsdf.anisotropic,vertex.uv,vertex.uv_screen_size,texture_pool);
    Real roughness = eval(bsdf.roughness,vertex.uv,vertex.uv_screen_size,texture_pool);
    Vector2  alpha = get_alpha(roughness,anisotropic);
//    alphaX = roughness * roughness * 2;
//    alphaY = roughness * roughness;


    //Geom item
    Real G = smith_masking_gtr2_anisotropy(to_local(frame, dir_in), alpha.x,alpha.y) *
            smith_masking_gtr2_anisotropy(to_local(frame, dir_out), alpha.x,alpha.y);
    Real D= get_D(to_local(frame,wh), alpha);

    //G = smith_masking_gtr2(to_local(frame, dir_in),roughness) * smith_masking_gtr2(to_local(frame, dir_out),roughness);
    return F * D *G / (4 * absDot(frame.n,dir_in));
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {


    if (dot(vertex.geometry_normal, dir_in) < 0 ||
            dot(vertex.geometry_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    Vector3  wh = normalize(dir_in+dir_out);
    // Homework 1: implement this!
    Real roughness = eval(bsdf.roughness,vertex.uv,vertex.uv_screen_size,texture_pool);
    Real anis = eval(bsdf.anisotropic,vertex.uv,vertex.uv_screen_size,texture_pool);
    auto alpha = get_alpha(roughness,anis);
    Real G= smith_masking_gtr2_anisotropy(to_local(frame,dir_in),alpha.x,alpha.y);
    return G * get_D(to_local(frame,wh),alpha) / (4 * dot(dir_in, frame.n));
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Vector3 local_dir_in = to_local(frame, dir_in);
    Real roughness = eval(
            bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real alpha = roughness * roughness;
    Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, alpha, rnd_param_uv);
    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{reflected,0.0,roughness};
    // Homework 1: implement this!

}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
