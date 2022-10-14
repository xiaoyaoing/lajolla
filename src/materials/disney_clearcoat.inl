#include "../microfacet.h"
#include "sstream"
static int count = 0;

Real get_alpha(Real gloss)
{
    return (1 - gloss) * 0.1 + gloss * 0.001;
}

Real get_Dc(Real alpha, Vector3 hl)
{
    Real alpha_sq = alpha * alpha;
    return (alpha_sq - 1) / (c_PI * log(alpha_sq)*(1 + (alpha_sq - 1)*(hl.z*hl.z)));
}


Real get_R0(Real eta)
{
    return pow(eta - 1,2) / pow(eta + 1,2);
}

Real get_G_clearcoat(Vector3 omega)
{
    double lambda = (sqrt(1 + (pow(omega.x * 0.25, 2) + pow(omega.y * 0.25, 2)) / (omega.z * omega.z)) - 1) / 2.0;
    return 1.0 / (1 + lambda);
}

Real get_Gc(Vector3 dir_in, Vector3 dir_out)
{
    return get_G_clearcoat(dir_in) * get_G_clearcoat(dir_out);
}

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    // Homework 1: implement this!

    Vector3 wh = normalize(dir_in+dir_out);
   // if( dot(wh,frame.n)<0) wh=-wh;

    //Compute Fresnel
    Real eta = 1.5;
    Real R0 = pow((eta-1)/(eta+1),2);
    Real f = R0 + (1-R0) * pow((1- absDot(wh,dir_out)),5);
    Spectrum Fc = Spectrum(f,f,f);
    //Compute D
    Real roughness = eval(bsdf.clearcoat_gloss,vertex.uv,vertex.uv_screen_size,texture_pool);
    Real alphaG = (1-roughness) * 0.1 + roughness * 0.001;
    Real D = (alphaG * alphaG - 1) / (M_PI * 2 * log(alphaG) *  (1+(alphaG * alphaG-1) * pow(to_local(frame,wh).z,2)) );
    Real G = smith_masking_gtr2(to_local(frame, dir_in), 0.25) *
             smith_masking_gtr2(to_local(frame, dir_out), 0.25);

//    if(rand()%10000 == 1)
//    {std::cout<<D<<" "<<G<<" "<<Fc<<" "<< to_local(frame,wh).z<<std::endl;
//      if(count++>50)
//        exit(-1);
//    }
    Real Dc = get_Dc(alphaG, to_local(frame, wh));
    Real F = get_R0(1.5) + (1 - get_R0(1.5)) * pow((1 - abs(dot(wh, dir_out))), 5);
    Real Gc = get_Gc(to_local(frame, dir_in), to_local(frame, dir_out));
    //std::cout << gloss << " " << alpha << " " << alpha << std::endl;
    //std::cout << Dc << " " << Fc << " " << Gc << std::endl;

  //  return make_zero_spectrum() + (F*Dc*Gc)/(4*abs(dot(frame.n,dir_in))x);
    return Fc * D * Gc / (4 * absDot(frame.n,dir_in));
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
  //  return max(dot(dir_out, frame.n), 0.0) / M_PI;

    Vector3 wh = normalize(dir_in + dir_out);
   // if( dot(wh,frame.n)<0) wh=-wh;

    auto cosTheta = to_local(frame,wh).z;
    auto cosTheta2 = cosTheta * cosTheta;

    Real roughness = eval(bsdf.clearcoat_gloss,vertex.uv,vertex.uv_screen_size,texture_pool);
    Real alphaG = (1-roughness) * 0.1 + roughness * 0.001;
    Real alphaG2 = alphaG * alphaG;


    auto res =(alphaG2-1)  / log(alphaG2)   / c_PI / (1+(alphaG2-1) * cosTheta2);
    if((res - get_Dc(alphaG, to_local(frame,wh)))>0.0001){
        auto t= get_Dc(alphaG, to_local(frame,wh));
        //std::cout<<t<<" "<<res<<"";
        throw("error");
                 get_Dc(alphaG, to_local(frame,wh));
    }



    return res * absDot(frame.n,wh) / (4 * absDot(wh,dir_out));
   // std::cout<<res<<" ";
  // return 1;
    return res;
    // Homework 1: implement this!

    return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }


    Real roughness = eval(bsdf.clearcoat_gloss,vertex.uv,vertex.uv_screen_size,texture_pool);
    Real alphaG = (1-roughness) * 0.1 + roughness * 0.001;
    assert(abs(alphaG-get_alpha(roughness))<0.0001);
    Real alphaG2 = alphaG * alphaG;
    Real cosTheta = sqrt( (1-pow(alphaG2,1-rnd_param_uv.x)) / (1- alphaG2) );
    Real sinTheta = sqrt(1-cosTheta * cosTheta);
    Real phi = rnd_param_uv.y * c_PI * 2;
    Vector3 wh = Vector3(sinTheta * cos(phi),sinTheta * sin(phi),cosTheta);

    Real  alpha =alphaG;
    Real h_ele = acos(sqrt((1 - pow((alpha * alpha), 1 - rnd_param_uv.x)) / (1 - alpha * alpha)));
    assert(abs(cos(h_ele)-cosTheta)<0.0001);


    wh = to_world(frame,wh);
    Vector3 out = normalize(-dir_in + 2*(dot(dir_in,wh)) * wh);
//    if(length(normalize((dir_in+out))-wh)>0.001){
//        std::cout<<dir_in<<" "<<out<<" "<<wh<<" "<<normalize(dir_in+out);
//        exit(-1);
//    }
    if(abs(length(out))-1>0.001){

    }
    auto rec = BSDFSampleRecord{
        out,
        //    to_world(frame,sample_cos_hemisphere(rnd_param_uv)),
        0,/*eta */
        alphaG
    };
   // std::cout<<rec.dir_out<<" "<<cosTheta<<" "<<sinTheta<<" "<<rnd_param_uv.x<<std::endl;
    if( isnan(rec.dir_out.x))
    {
        std::stringstream  s;
        s<<dir_in<<" "<<wh<<" "<<rec.dir_out<<" "<<cosTheta<<" "<<frame.n;
        std::cout<<s.str();
        assert(false);
    }
    return rec;
    // Homework 1: implement this!

    return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
