#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    if(vertex.explictRefract )
        reflect = !vertex.refract;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum baseColor = eval(bsdf.base_color,vertex.uv,vertex.uv_screen_size,texture_pool);
    Real roughness  = eval(bsdf.roughness,vertex.uv,vertex.uv_screen_size,texture_pool);
    //roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real eta = dot(vertex.geometry_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Vector3 wh;
    if (reflect) {
        wh = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        wh = normalize(dir_in + dir_out * eta);
    }

    if (dot(wh, frame.n) < 0) {
        wh = -wh;
    }

    Real  Fg = fresnel_dielectric(dot(wh,dir_in),eta);

    Real  D = GTR2(dot(wh,frame.n),roughness);
    frame = Frame(wh);
    Real  G = smith_masking_gtr2(to_local(frame, dir_in),roughness) *
              smith_masking_gtr2(to_local(frame, dir_out),roughness);
    // Reflect case
    if(reflect){
        //assert(lum<1000) ;
        return baseColor * Fg * D * G  / (4 *
                absDot(dir_in,frame.n));
    }
    // Refract Case
    else {
        Real eta_factor = dir == TransportDirection::TO_LIGHT ? (1 / (eta * eta)) : 1;
        //eta_factor = 1;
        Real h_dot_out = dot(wh,dir_out);
        Real h_dot_in = dot(wh,dir_in);
        Real deom = h_dot_in+eta * h_dot_out;
        Real F = Fg;
        Spectrum res = baseColor * (eta_factor * (1 - F) * D * G  * fabs(h_dot_out * h_dot_in)) /
        (fabs(dot(frame.n, dir_in)) * deom * deom );
        if(abs(deom)<0.001)
            std::cout<<deom<<std::endl;
        {

//            std::cout<<D<<" "<<F<<" "<<G<<" "<<eta_factor * eta<<" "<< abs(h_dot_out * h_dot_in) /
//            fabs(dot(frame.n, dir_in))<<std::endl;
        }
        return  sqrt(baseColor) * (1-Fg) * D * G * abs(h_dot_out * h_dot_in) /
                 ( absDot(frame.n,dir_in) * abs(deom * deom) ) ;
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    if(vertex.explictRefract )
        reflect = !vertex.refract;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real eta = dot(vertex.geometry_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real roughness = eval(
            bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);
    Real D = GTR2(dot(half_vector, frame.n), roughness);
    Real G_in = smith_masking_gtr2(to_local(frame, dir_in), roughness);
    if (reflect) {
       // std::cout<<(F * D * G_in) / (4 * fabs(dot(frame.n, dir_in)))<<" ";
        return (F * D * G_in) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        Real h_dot_out = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);

     //   std::cout<< (1 - F) * D * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in))<<" ";

        return (1 - F) * D * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    Real eta = dot(vertex.geometry_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real roughness = eval(
            bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Sample a micro normal and transform it to world space -- this is our half-vector.
    Real alpha = roughness * roughness;
    Vector3 local_dir_in = to_local(frame, dir_in);
    Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(half_vector, dir_in);

//    Real  Rs = (h_dot_in - eta * h_dot_out_temp) /(h_dot_in + eta * h_dot_out_temp);
//    Real  Rp = (eta * h_dot_in -  h_dot_out_temp) / (eta * h_dot_in + h_dot_out_temp);
//    Real  Fg = 0.5 * (Rs * Rs) / (Rp * Rp);
    Real F = fresnel_dielectric(h_dot_in, eta);

    if (rnd_param_w <= F) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness,false};
    } else {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            return {};
        }
        // flip half_vector if needed
        if (h_dot_in < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_out= sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
        assert(dot(dir_in,half_vector) * dot(refracted,half_vector)<0);
        return BSDFSampleRecord{refracted, eta, roughness,true};
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
