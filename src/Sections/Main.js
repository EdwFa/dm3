import React, { Component } from "react";
import { useState, useEffect, createRef } from "react";

import { Navigate, Link } from "react-router-dom";

import Graph from "react-graph-vis";
import { v4 as uuidv4 } from "uuid";
//import "./network.css";

import { AgGridReact } from "ag-grid-react";
import "ag-grid-enterprise";
import "ag-grid-community/styles/ag-grid.css";
import "ag-grid-community/styles/ag-theme-alpine.css";

import { variables } from "./Variables.js";

export class Main extends Component {
  constructor(props) {
    super(props);

    this.gridRef = createRef();
    this.gridAnaliseRef = createRef();
    this.state = {
      token: variables.token,
    };
  }

  componentDidMount() {
    console.log("start");
  }

  render() {
    const { token } = this.state;

    if (!token) {
      return <Navigate push to="/login" />;
    } else {
      return (
        <>
          <div className="bg-zinc-100 py-20 h-screen">
            <div className="mx-24">
              <div className="flex justify-center">
                <div className="grid grid-cols-2 gap-4">
                  <div>
                    <div class="max-w-sm p-6 bg-white border border-gray-200 rounded-lg shadow mb-4">
                      <Link to="/tematic_review">
                        <h5 class="mb-2 text-2xl font-bold tracking-tight text-gray-900">
                          Тематическое моделирование
                        </h5>
                      </Link>
                      <p class="mb-3 font-normal text-gray-700">
                        Подсистема для извлечения фактов знаний, оснащенная
                        фильтро-поисковой системой.
                      </p>
                      <a
                        href="/tematic_review"
                        class="inline-flex items-center px-3 py-2 text-sm font-medium text-center text-white bg-blue-700 rounded-lg hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300"
                      >
                        Начать поиск
                        <svg
                          class="w-3.5 h-3.5 ml-2"
                          aria-hidden="true"
                          xmlns="http://www.w3.org/2000/svg"
                          fill="none"
                          viewBox="0 0 14 10"
                        >
                          <path
                            stroke="currentColor"
                            stroke-linecap="round"
                            stroke-linejoin="round"
                            stroke-width="2"
                            d="M1 5h12m0 0L9 1m4 4L9 9"
                          />
                        </svg>
                      </a>
                    </div>
                    <div class="max-w-sm p-6 bg-white border border-gray-200 rounded-lg shadow">
                      <Link to="/ddi_review">
                        <h5 class="mb-2 text-2xl font-bold tracking-tight text-gray-900">
                          Поиск ключевых слов
                        </h5>
                      </Link>
                      <p class="mb-3 font-normal text-gray-700">
                        Подсистема для извлечения фактов знаний, оснащенная
                        фильтро-поисковой системой.
                      </p>
                      <a
                        href="/ddi_review"
                        class="inline-flex items-center px-3 py-2 text-sm font-medium text-center text-white bg-blue-700 rounded-lg hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 dark:bg-blue-600"
                      >
                        Начать поиск
                        <svg
                          class="w-3.5 h-3.5 ml-2"
                          aria-hidden="true"
                          xmlns="http://www.w3.org/2000/svg"
                          fill="none"
                          viewBox="0 0 14 10"
                        >
                          <path
                            stroke="currentColor"
                            stroke-linecap="round"
                            stroke-linejoin="round"
                            stroke-width="2"
                            d="M1 5h12m0 0L9 1m4 4L9 9"
                          />
                        </svg>
                      </a>
                    </div>
                  </div>
                  <div className="bg-gray-300"></div>
                </div>
              </div>
            </div>
          </div>
        </>
      );
    }
  }
}