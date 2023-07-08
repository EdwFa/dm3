import React, { Component } from 'react';
import { useState, useEffect, createRef } from 'react';

import { Navigate, Link } from 'react-router-dom';

import Graph from "react-graph-vis";
import { v4 as uuidv4 } from 'uuid'
//import "./network.css";

import { AgGridReact } from 'ag-grid-react';
import 'ag-grid-enterprise';
import 'ag-grid-community/styles/ag-grid.css';
import 'ag-grid-community/styles/ag-theme-alpine.css';

import { variables } from './Variables.js';


export class Main extends Component {

    constructor(props) {
        super(props);

        this.gridRef = createRef();
        this.gridAnaliseRef = createRef();
        this.state = {
            token: variables.token,
        }
    }

    componentDidMount() {
        console.log('start');
    }

    render() {
        const {
            token,
        } = this.state;

        if (!token){
            return <Navigate push to="/login" />
        } else {
            return (
            <>
                <ul class="pt-4 mt-4 space-y-2 border-t border-gray-200 dark:border-gray-700">
                    <Link to="/tematic_review">
                      <li>
                        <a href="#" class="flex items-center p-2 text-base font-normal text-gray-900 rounded-lg dark:text-white hover:bg-gray-100 dark:hover:bg-gray-700">
                          <span class="flex-1 ml-3 whitespace-nowrap">Тематическое моделирование</span>
                        </a>
                      </li>
                    </Link>
                    <Link to="/ddi_review">
                      <li>
                        <a href="#" class="flex items-center p-2 text-base font-normal text-gray-900 rounded-lg dark:text-white hover:bg-gray-100 dark:hover:bg-gray-700">
                          <span class="flex-1 ml-3 whitespace-nowrap">Поиск ключевых слов</span>
                        </a>
                      </li>
                    </Link>
                </ul>
            </>
            )
        }
    }
}