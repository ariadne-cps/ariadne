/***************************************************************************
 *            concurrency/concurrency_manager.cpp
 *
 *  Copyright  2007-20  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <time.h>

#include "../utility/macros.hpp"
#include "../concurrency/concurrency_manager.hpp"

namespace Ariadne {

ConcurrencyManager::ConcurrencyManager() : _maximum_concurrency(std::thread::hardware_concurrency()), _concurrency(1) {}

SizeType ConcurrencyManager::maximum_concurrency() const {
    return _maximum_concurrency;
}

SizeType ConcurrencyManager::concurrency() const {
    return _concurrency;
}

void ConcurrencyManager::set_concurrency(SizeType value) {
    ARIADNE_PRECONDITION(value <= _maximum_concurrency and value > 0);
    std::lock_guard<std::mutex> lock(_data_mutex);
    _concurrency = value;
}

List<TaskSearchPointAppraisal> ConcurrencyManager::last_search_best_points() const {
    return _last_search_best_points;
}

void ConcurrencyManager::set_last_search_best_points(List<TaskSearchPointAppraisal> const &points) {
    std::lock_guard<std::mutex> lock(_data_mutex);
    _last_search_best_points = points;
}

auto ConcurrencyManager::last_property_refinement_values() const -> PropertyRefinementsMap {
    return _last_property_refinement_values;
}

void ConcurrencyManager::set_last_property_refinement_values(PropertyRefinementsMap const &refinements) {
    _last_property_refinement_values.clear();
    _last_property_refinement_values.adjoin(refinements);
}

void ConcurrencyManager::print_last_search_best_points() const {
    if (not _last_search_best_points.empty()) {
        std::ofstream file;
        file.open("points.m");
        auto space = _last_search_best_points.front().point().space();
        auto dimension = space.dimension();
        auto size = _last_search_best_points.size();
        file << "x = [1:" << size << "];\n";
        Map<SizeType,List<int>> values;
        for (SizeType i=0; i<dimension; ++i) values.insert(make_pair(i,List<int>()));
        for (auto appraisal : _last_search_best_points) {
            auto point = appraisal.point();
            for (SizeType i=0; i<dimension; ++i) values.at(i).push_back(point.coordinates()[i]);
        }
        file << "figure(1);\n";
        file << "hold on;\n";
        for (SizeType i=0; i<dimension; ++i) {
            file << "y" << i << " = " << values.at(i) << ";\n";
            auto name = space.parameters()[i].path().last();
            std::replace(name.begin(), name.end(), '_', ' ');
            file << "plot(x,y" << i << ",'DisplayName','" << name << "');\n";
        }
        file << "legend;\n";
        file << "hold off;\n";
        file.close();
    }
}

void ConcurrencyManager::print_last_property_refinement_values() const {
    if (not _last_property_refinement_values.empty()) {
        for (auto prop : _last_property_refinement_values) {
            std::ofstream file;
            file.open("refinements.m");
            List<ExactDouble> values;
            for (auto p : prop.second) {
                values.push_back(dynamic_cast<RangeConfigurationProperty<ExactDouble>*>(p.get())->get());
            }
            file << "y = " << values << ";\n";
            file << "x = [1:" << values.size() << "];\n";
            file << "semilogy(x,y)";
            file.close();
        }
    }
}

}

inline bool _init_randomiser() {
    srand(time(nullptr));
    return true;
}

static const bool init_randomiser = _init_randomiser();