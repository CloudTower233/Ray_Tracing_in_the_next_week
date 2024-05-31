#pragma once
#ifndef CONSTANT_MEDIUM_H
#define CONSTANT_MEDIUM_H

#include "rtweekend.h"

#include "hittable.h"
#include "material.h"
#include "texture.h"

class constant_medium : public hittable {
public:
    constant_medium(shared_ptr<hittable> b, double d, shared_ptr<texture> a)
        : boundary(b),
        neg_inv_density(-1 / d),
        phase_function(make_shared<isotropic>(a))
    {}

    constant_medium(shared_ptr<hittable> b, double d, color c)
        : boundary(b),
        neg_inv_density(-1 / d),
        phase_function(make_shared<isotropic>(c))
    {}

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        return boundary->bounding_box(time0, time1, output_box);
    }

public:
    shared_ptr<hittable> boundary;        // �����İ�Χ��
    shared_ptr<material> phase_function;  // ����
    double neg_inv_density;  // Ϊ�˷�����㣬�ܶ�ʹ�õ��������ڸ��ʺ������õ�log����������Ҳ���ϸ���
};

bool constant_medium::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    // Print occasional samples when debugging. To enable, set enableDebug true.
    const bool enableDebug = false;
    const bool debugging = enableDebug && random_double() < 0.00001;

    hit_record rec1, rec2;

    // �����ж�����������˶��ӳ������޳����Ƿ�������������İ�Χ���ཻ
    if (!boundary->hit(r, -infinity, infinity, rec1))
        // ������ཻ�򷵻�false
        return false;

    // �������ڵ�һ���ཻ�� ��һ���ӳ������޳����Ƿ��������ڶ����ཻ��������ܵڶ����ཻ�򷵻�false
    if (!boundary->hit(r, rec1.t + 0.0001, infinity, rec2))
        return false;

    if (debugging) std::cerr << "\nt_min=" << rec1.t << ", t_max=" << rec2.t << '\n';

    // �������ཻ��ʱ�������� t_min��t_max֮��
    if (rec1.t < t_min) rec1.t = t_min;
    if (rec2.t > t_max) rec2.t = t_max;

    // �����һ���ཻ��ʱ��ȵڶ��λ�Ҫ���򷵻�false
    if (rec1.t >= rec2.t)
        return false;


    // �ѵ�һ���ཻ��ʱ��������0���ϣ������ߵ������
    if (rec1.t < 0)
        rec1.t = 0;

    // �����ܳ���
    const auto ray_length = r.direction().length();
    // ������������ڵĳ���
    const auto distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
    // ���һ�������� ������ �ܶ� ��������һ����ײ�����
    const auto hit_distance = neg_inv_density * log(random_double());

    // �����ǰ������������ڵĳ��� �� ��ǰ����������ײ����� ��Ҫ�̣��򷵻�false
    if (hit_distance > distance_inside_boundary)
        return false;

    // ��ײʱ�� = ��һ����ײ���ʱ��t + �����/�����ܳ���
    rec.t = rec1.t + hit_distance / ray_length;
    rec.p = r.at(rec.t);

    if (debugging) {
        std::cerr << "hit_distance = " << hit_distance << '\n'
            << "rec.t = " << rec.t << '\n'
            << "rec.p = " << rec.p << '\n';
    }

    rec.normal = vec3(1, 0, 0);  // ��ײ��������
    rec.front_face = true;     // ����ײ��ǰ����
    rec.mat_ptr = phase_function;

    return true;
}
#endif